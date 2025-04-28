#' @title Replace internal standards in SummarizedExperiment
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param internalStandards
#' @importFrom dplyr %>%
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
replaceInternalStandards <- function(exp, internalStandards = rowData(exp)$suggestedIS) {
    if (is.null(internalStandards)) {
        return(exp)
    }
    message("Replacing Internal Standards")
    is <- rowData(exp)$compound_is
    indexes <- match(internalStandards, is)

    toReplace <- which(is.na(indexes))
    indexes[toReplace] <- toReplace

    rowData(exp)$compound_is <- is[indexes]
    area_is <- metadata(exp)$secondary

    assay(exp, area_is, withDimnames = FALSE) <- assay(exp, area_is)[indexes, ]

    # exp <- calculateRatio(exp) %>%
    #     addBatchCorrection()

    return(exp)
}

#' @title Calculate the RSDQC of an assay
#' @description The RSDQC is an indicator of variance of compounds in QC
#' samples. This function returns the RSDQC of all compounds in the
#' SummarizedExperiment
#' @details placeholder
#' @param exp SummarizedExperiment object
#' @param assay Assay to use for calculating the RSDQC. Defaults to
#' "ratio_corrected"
#' @param type
#' @importFrom matrixStats rowSds
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' rsdqc(exp)
rsdqc <- function(exp, assay = "ratio_corrected", type = metadata(exp)$QC) {
    qc_df <- assay(exp[, exp$type == type], assay)
    sds <- rowSds(qc_df, na.rm = TRUE)
    means <- rowMeans(qc_df, na.rm = TRUE)

    return(sds / means * 100)
}


#' @title Calculate the RSDQC of an assay
#' @description The RSDQC is an indicator of variance of compounds in QC
#' samples. This function returns the RSDQC of all compounds in the
#' SummarizedExperiment
#' @details placeholder
#' @param exp SummarizedExperiment object
#' @importFrom matrixStats rowSds
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' calculateCorrectedRSDQCs2(exp)
calculateCorrectedRSDQCs2 <- function(exp, primaryAssay = "area",
                                      secondaryAssay = "area_is",
                                      qcType = "SQC", returnRSDs = FALSE,
                                      saveAssay = "ratio_corrected") {
    if (metadata(exp)$hasIS) {
        comp_is <- unique(rowData(exp)$compound_is)
        df <- expand.grid(
            compound_is = comp_is,
            compound = rownames(exp)
        )

        comp_row <- match(df$compound_is, rowData(exp)$compound_is)

        areas <- assay(exp, primaryAssay)[df$compound, ]
        area_is <- assay(exp, secondaryAssay)[comp_row, ]
        ratio <- areas / area_is
    } else {
        comp_is <- seq_len(nrow(exp))
        ratio <- assay(exp, primaryAssay)
    }




    dims <- dimnames(ratio)

    # Within Batch correction
    withinRatio <- do.call(cbind, lapply(unique(exp$batch), function(j) {
        idx <- exp$batch == j
        subExp <- exp[, idx]

        x <- ratio[, idx]

        order <- which(subExp$type == metadata(exp)$QC)
        vals <- x[, order]

        orderMatrix <- matrix(
            data = rep(order, nrow(vals)),
            ncol = length(order),
            byrow = TRUE
        )
        slope <- rowWiseSlope(orderMatrix, vals)
        int <- rowWiseIntercept(orderMatrix, y = vals, slope = slope)
        orderMatrix <- matrix(
            data = rep(seq_len(ncol(x)), nrow(vals)),
            ncol = ncol(x),
            byrow = TRUE
        )
        vals <- orderMatrix * slope + int
        vals
    }))


    m <- ratio / withinRatio

    m <- do.call(cbind, lapply(unique(exp$batch), function(j) {
        idx <- exp$batch == j
        m[, idx] * rowMedians(ratio[, idx], na.rm = TRUE)
    }))


    m[!is.finite(m)] <- NA
    m[m < 0] <- NA


    dimnames(m) <- dims
    ratio <- m

    rm(m)



    # subset the experiment with only QCs
    qcExp <- exp[, exp$type == qcType]


    # Between Batch -----------------------------------------------------------


    compound_qc_ratio_median <- rowMedians(ratio[, colnames(qcExp)], na.rm = TRUE)

    # For each batch, calculate the correction factor and multiply the
    # Ratio to obtain the Ratio Corrected factors
    ratio <- do.call(cbind, lapply(unique(exp$batch), function(batch) {
        correction_factor <- 1

        # Check if the batch is present in the QC-subsetted experiment
        if (batch %in% qcExp$batch) {
            # Subset the QC experiment for the current batch
            qcBatchExp <- qcExp[, qcExp$batch == batch]


            # Calculate the median Ratio for this batch
            med_ratio_qc <- rowMedians(ratio[, colnames(qcBatchExp)], na.rm = TRUE)

            # Calculate the QC correction factor by dividing the overall median
            # by this batch median
            correction_factor <- compound_qc_ratio_median / med_ratio_qc
        }
        # Multiply the entire batch with the calculated correction factors
        ratio[, colnames(exp[, exp$batch == batch])] * correction_factor
    }))[, colnames(exp)]


    # Ensure that the right column names are used, in the right order
    qc_ratio <- ratio[, colnames(exp[, exp$type == qcType]), drop = FALSE]
    sds <- rowSds(qc_ratio, na.rm = TRUE)
    means <- rowMeans(qc_ratio, na.rm = TRUE)
    rsd <- sds / means * 100


    res <- split(rsd, ceiling(seq_along(rsd) / length(comp_is)))

    m <- do.call(cbind, res)
    dimnames(m) <- list(comp_is, rownames(exp))
    return(m)
}
