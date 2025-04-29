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
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' matrixRSDQCs(exp)
matrixRSDQCs <- function(exp,
                         primaryAssay = metadata(exp)$primary,
                         secondaryAssay = metadata(exp)$secondary,
                         qcType = metadata(exp)$QC) {

    if (metadata(exp)$hasIS) {
        comp_is <- unique(rowData(exp)$compound_is)
        df <- expand.grid(
            compound_is = comp_is,
            compound = rownames(exp)
        )

        comp_row <- match(df$compound_is, rowData(exp)$compound_is)

        areas <- assay(exp, primaryAssay)[df$compound, ]
        area_is <- assay(exp, secondaryAssay)[comp_row, ]
        ratios <- areas / area_is
    } else {
        comp_is <- seq_len(nrow(exp))
        ratios <- assay(exp, primaryAssay)
    }

    qcs <- exp$type == qcType
    ratios <- betweenBatchCorrection(
        ratio = ratios,
        qcAliquots = colnames(exp)[qcs],
        qcBatches = exp$batch[qcs],
        aliquots = colnames(exp),
        batches = exp$batch
    )

    # Ensure that the right column names are used, in the right order
    qc_ratio <- ratios[, colnames(exp[, exp$type == qcType]), drop = FALSE]
    sds <- rowSds(qc_ratio, na.rm = TRUE)
    means <- rowMeans(qc_ratio, na.rm = TRUE)
    rsd <- sds / means * 100

    res <- split(rsd, ceiling(seq_along(rsd) / length(comp_is)))
    m <- do.call(cbind, res)
    dimnames(m) <- list(comp_is, rownames(exp))
    return(m)
}
