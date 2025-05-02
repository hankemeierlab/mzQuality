#' @title Replace internal standards in SummarizedExperiment
#' @description Replaces the internal standards in a SummarizedExperiment
#'     object with the specified internal standards or the suggested ones.
#' @details This function updates the internal standards (`compound_is`) in
#'     the `rowData` of the SummarizedExperiment object. If no internal
#'     standards are provided, the function uses the suggested internal
#'     standards (`rowData(exp)$suggestedIS`). It also updates the secondary
#'     assay (`area_is`) to reflect the new internal standards.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param internalStandards A character vector specifying the new internal
#'     standards to use. Defaults to the suggested internal standards
#'     (`rowData(exp)$suggestedIS`).
#' @return A SummarizedExperiment object with updated internal standards and
#'     secondary assay values.
#' @importFrom dplyr %>%
#' @examples
#' # Example usage:
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Replace internal standards with suggested ones
#' exp <- replaceInternalStandards(exp)
#'
#' # Replace internal standards with a randomized vector of already present
#' # internal standards
#' customIS <- sample(rowData(exp)$compound_is)
#' exp <- replaceInternalStandards(exp, internalStandards = customIS)
#' @export
replaceInternalStandards <- function(
        exp, internalStandards = rowData(exp)$suggestedIS
) {
    if (is.null(internalStandards)) {
        return(exp)
    }

    is <- rowData(exp)$compound_is
    indexes <- match(internalStandards, is)

    toReplace <- which(is.na(indexes))
    indexes[toReplace] <- toReplace

    rowData(exp)$compound_is <- is[indexes]
    area_is <- metadata(exp)$secondary

    assay(exp, area_is, withDimnames = FALSE) <- assay(exp, area_is)[indexes, ]

    return(exp)
}

#' @title Calculate the RSDQC of an assay
#' @description Computes the Relative Standard Deviation for Quality Control
#'     (RSDQC) of compounds in QC samples from a SummarizedExperiment object.
#' @details The RSDQC is a measure of the variance of compounds in QC samples,
#'     expressed as a percentage. This function calculates the RSDQC for all
#'     compounds in the specified assay of the SummarizedExperiment object.
#'     It uses the standard deviation and mean of the QC sample values to
#'     compute the RSDQC. Missing values are ignored during the calculation.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param assay A string specifying the assay to use for calculating the
#'     RSDQC. Defaults to `"ratio_corrected"`.
#' @param type A string specifying the type of QC samples to use for
#'     calculating the RSDQC. Defaults to the `QC` metadata of the experiment.
#' @return A numeric vector containing the RSDQC values for all compounds in
#'     the specified assay.
#' @importFrom matrixStats rowSds
#' @examples
#' # Example usage:
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Calculate the RSDQC for the default assay and QC type
#' rsdqc_values <- rsdqc(exp)
#'
#' # Calculate the RSDQC for a specific assay and QC type
#' rsdqc_values <- rsdqc(exp, assay = "ratio", type = "QC")
#' @export
rsdqc <- function(exp, assay = "ratio_corrected", type = metadata(exp)$QC) {
    qc_df <- assay(exp[, exp$type == type], assay)
    sds <- rowSds(qc_df, na.rm = TRUE)
    means <- rowMeans(qc_df, na.rm = TRUE)
    return(sds / means * 100)
}


#' @title Calculate the RSDQC of an assay
#' @description Computes the Relative Standard Deviation for Quality Control
#'     (RSDQC) of compounds in QC samples from a SummarizedExperiment object.
#' @details This function calculates the RSDQC, which is an indicator of the
#'     variance of compounds in QC samples. It supports both primary and
#'     secondary assays and performs between-batch correction if internal
#'     standards (IS) are available. The result is a matrix of RSDQC values
#'     for each compound and internal standard.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param primaryAssay A string specifying the name of the primary assay to
#'     use. Defaults to the `primary` metadata of the experiment.
#' @param secondaryAssay A string specifying the name of the secondary assay
#'     to use. Defaults to the `secondary` metadata of the experiment.
#' @param qcType A string specifying the type of QC samples. Defaults to the
#'     `QC` metadata of the experiment.
#' @return A matrix of RSDQC values, with rows corresponding to internal
#'     standards and columns corresponding to compounds.
#' @importFrom matrixStats rowSds
#' @examples
#' # Example usage:
#' # exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#' # rsdMatrix <- matrixRSDQCs(exp)
#' @export
#' @examples
#' # Example usage:
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Calculate the RSDQC matrix for the default primary and secondary assays
#' rsd_matrix <- matrixRSDQCs(exp)
matrixRSDQCs <- function(
        exp, primaryAssay = metadata(exp)$primary,
        secondaryAssay = metadata(exp)$secondary,
        qcType = metadata(exp)$QC
) {

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
    ratios <- .betweenBatchCorrection(
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
