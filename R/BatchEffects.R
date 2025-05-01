#' @title Between-Batch Correction for Metabolomics Data
#' @description Normalizes ratio data across multiple batches using QC samples
#'   as reference points.
#' @details This function calculates batch-specific correction factors based on
#'   QC sample measurements. For each batch, it computes the median ratio of QC
#'   samples within the batch and divides the global QC median by this
#'   batch-specific median to obtain a correction factor. The correction factor
#'   is then applied to all samples within that batch to normalize their values
#'   relative to other batches.
#'
#'   If a batch contains no QC samples, a correction factor of 1 is used
#'   (no correction).
#'
#' @param ratio Matrix of ratio values with compounds as rows and samples as
#' columns
#' @param qcAliquots Character vector of column names identifying QC samples
#' @param qcBatches Vector of batch identifiers for the QC samples
#' @param aliquots Character vector of all sample identifiers (column names)
#' @param batches Vector of batch identifiers for all samples
#'
#' @return A matrix of the same dimensions as the input `ratio`, containing
#'   the batch-corrected values
#'
#' @importFrom matrixStats rowMedians
#' @keywords internal
betweenBatchCorrection <- function(
        ratio, qcAliquots, qcBatches, aliquots, batches
) {

    compound_qc_ratio_median <- rowMedians(ratio[, qcAliquots], na.rm = TRUE)

    # For each batch, calculate the correction factor and multiply the
    # Ratio to obtain the Ratio Corrected factors
    for (batch in unique(batches)) {
        correction_factor <- 1

        # Check if the batch is present in the QC-subsetted experiment
        if (batch %in% qcBatches) {
            # Subset the QC experiment for the current batch
            qcAliqs <- qcAliquots[qcBatches == batch]

            # Calculate the median Ratio for this batch
            med_ratio_qc <- rowMedians(ratio[, qcAliqs], na.rm = TRUE)

            # Calculate the QC correction factor by dividing the overall median
            # by this batch median
            correction_factor <- compound_qc_ratio_median / med_ratio_qc

            # Multiply the entire batch with the calculated correction factors
            aliqs <- aliquots[batches == batch]
            ratio[, aliqs] <- ratio[, aliqs] * correction_factor
        }
    }

    return(ratio)
}

#' @title Batch Effect Correction
#' @description Applies batch-wise correction to normalize ratios across
#'   multiple batches using QC samples.
#' @details This function uses QC sample medians per batch and across all
#'   batches to calculate correction factors. These are used to normalize
#'   the specified assay across all batches, adjusting for batch effects.
#'   Optionally, outliers identified via `exp$use` can be removed before
#'   calculating corrections.
#' @param exp A `SummarizedExperiment` object with metabolomics data,
#'   including batch and type metadata.
#' @param assayName Character string indicating which assay contains the
#'   raw ratio values (default is "ratio").
#' @param saveAssay Character string indicating where to store corrected
#'   results (default is "ratio_corrected").
#' @param qc Character string denoting the QC sample type. Defaults to
#'   `metadata(exp)$QC`.
#' @param removeOutliers Logical flag to exclude outlier samples using the
#'   `use` column in `colData(exp)` (default is TRUE).
#' @returns A `SummarizedExperiment` object with batch-corrected ratios in
#'   a new assay slot.
#' @importFrom matrixStats rowMedians
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Add batch correction assays
#' exp <- addBatchCorrectionAssay(exp)
#'
#' # Show results of the assay
#' assay(exp, "ratio_corrected")
addBatchCorrectionAssay <- function(
        exp, assayName = "ratio", saveAssay = "ratio_corrected",
        qc = metadata(exp)$QC, removeOutliers = TRUE
) {
    # check if the experiment is valid
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    ratios <- assay(exp, assayName)
    if (removeOutliers) {
        ratios[, !exp$use] <- NA
    }

    qcs <- exp$type == qc
    ratios <- betweenBatchCorrection(
        ratio = ratios,
        qcAliquots = colnames(exp)[qcs],
        qcBatches = exp$batch[qcs],
        aliquots = colnames(exp),
        batches = exp$batch
    )

    ratios[!is.finite(ratios)] <- NA
    ratios[, !exp$use] <- NA
    assay(exp, saveAssay) <- ratios

    # Return the experiment
    return(exp)
}

#' @title Perform batch-wise rotation of assay values
#' @description Adjusts assay values within a batch using a linear model
#'     constructed from quality control (QC) samples.
#' @details This function performs a batch-wise rotation of assay values to
#'     account for systematic variations within a batch. It uses QC samples to
#'     construct a linear model (slope and intercept) and applies this model to
#'     adjust the assay values for all samples in the batch. Values with too
#'     many missing data points are excluded from rotation.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param orders A numeric vector specifying the injection order of the samples.
#' @param assay A numeric matrix or data frame containing the assay values to
#'     be adjusted.
#' @param qcColumns A logical vector indicating which columns correspond to QC
#'     samples.
#' @param batch A character or numeric value specifying the batch to process.
#' @param withinBatchNonNA An integer specifying the minimum number of non-NA
#'     values required within a batch for reliable rotation. Defaults to 3.
#' @return A numeric matrix with the rotated assay values, adjusted for
#'     systematic batch effects.
#' @importFrom stats median
#' @examples
#' # Example usage:
#' # exp <- SummarizedExperiment(...)
#' # orders <- c(1, 2, 3, 4, 5)
#' # assay <- matrix(rnorm(50), nrow = 10)
#' # qcColumns <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
#' # batch <- "Batch1"
#' # rotatedAssay <- performRotation(exp, orders, assay, qcColumns, batch)
performRotation <- function(
        exp, orders, assay, qcColumns,
        batch, withinBatchNonNA = 3
){

    N <- nrow(exp)
    # retrieve the data from the current batch
    batchColumns <- exp$batch == batch
    qcBatchColumns <- batchColumns & qcColumns

    # construct order matrix using qc samples
    qcBatchOrders <- orders[qcBatchColumns]


    orderMatrix <- matrix(
        data = rep(qcBatchOrders, N),
        ncol = sum(qcBatchColumns),
        byrow = TRUE
    )

    # Calculate the slope and intercept for the linear model
    qcVals <- assay[, qcBatchColumns, drop = FALSE]
    slope <- rowWiseSlope(orderMatrix, qcVals)
    int <- rowWiseIntercept(orderMatrix, y = qcVals, slope = slope)

    orderMatrix <- matrix(
        data = rep(orders[batchColumns], N),
        ncol = sum(batchColumns),
        byrow = TRUE
    )

    # Calculate the rotation matrix
    vals <- assay[, batchColumns, drop = FALSE]
    rotation <- orderMatrix * slope + int
    batchMedians <- apply(vals, 1, median, na.rm = TRUE)

    # Fix values with too many NAs (unreliable model, should not rotate)
    tooManyNAs <- rowSums(!is.na(vals)) < withinBatchNonNA

    rotation[tooManyNAs, ] <- 1
    batchMedians[tooManyNAs] <- 1

    return((vals / rotation) * batchMedians)
}


#' @title Within-Batch Signal Drift Correction
#' @description Applies within-batch signal correction to assay data based
#'   on QC sample drift using a linear model per feature.
#' @details This function corrects for within-batch signal drift in
#'   metabolomics data by fitting a linear model to quality control (QC)
#'   samples across the injection order within each batch. Outliers in the
#'   QC measurements at the beginning and end of each batch can optionally
#'   be removed before model fitting. The correction is applied to the
#'   target assay by dividing the original intensities by the predicted
#'   values and rescaling using the median QC intensities.
#' @returns A `SummarizedExperiment` object with a new assay containing
#'   the drift-corrected values.
#'
#' @param exp A `SummarizedExperiment` object containing metabolomics data
#'   with `colData` columns: `batch`, `order`, `type`, and `injection_time`.
#' @param assayName Character string specifying the name of the assay to
#' correct (default is "ratio").
#' @param saveAssay Character string specifying the name of the new assay to
#'   store corrected values (default is "ratio_corrected").
#' @param qcType Character string specifying the sample type used for QC
#'   (usually "SQC" or "QC"). Defaults to the value stored in
#'   `metadata(exp)$QC`.
#' @param withinBatchNonNA Integer specifying the minimum number of
#' non-NA values before a correction is applied (default is 3).
#'
#' @importFrom dplyr group_by row_number pull n mutate arrange filter
addWithinBatchCorrectionAssay <- function(
        exp, assayName = "ratio", saveAssay = "ratio_corrected",
        qcType = metadata(exp)$QC, withinBatchNonNA = 3
) {

    # Multiplying factor which prevents negative values when the slope
    # goes through zero
    factor <- 1e4

    N <- nrow(exp)
    assay <- assay(exp, assayName)

    # set all batches
    batches <- unique(exp$batch)

    orders <- do.call(c, lapply(batches, function(batch) {
        batchColumns <- exp$batch == batch
        if ("injection_time" %in% colnames(colData(exp))) {
            orders <- order(exp$injection_time[batchColumns])
        } else {
            orders <- seq_along(which(batchColumns))
        }
        orders
    }))


    qcColumns <- exp$type == qcType

    for (batch in batches) {
        rotated <- performRotation(
            exp = exp,
            orders = orders,
            assay = assay * factor,
            qcColumns = qcColumns,
            batch = batch,
            withinBatchNonNA = withinBatchNonNA
        )
        rotated[!is.finite(rotated)] <- NA
        assay[, exp$batch == batch] <- rotated / factor
    }

    assay(exp, saveAssay) <- assay
    return(exp)
}


#' @title Add batch correction assays and metrics
#' @description Batch effects can cause the variation of A/IS ratios to be
#' very high. By using pooled QC samples, batch effect factors per compound
#' can be calculated. These are multiplied by the original A/IS ratio to obtain
#' batch-corrected ratios. This usually results in improved results.
#' @param exp SummarizedExperiment Object
#' @param assay Assay to be used for batch correction
#' @param qcType Type of QC to be used for batch correction
#' @param removeOutliers Remove outliers from the batch correction
#' @param useWithinBatch Use within batch correction
#' @param calculateRSDQC Boolean value. Should the assay be used to
#' calculate RSDQC values? Default is FALSE.
#' @importFrom matrixStats rowSds
#' @returns SummarizedExperiment with batch-effect corrected RSDQCs in the
#' rowData slot
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Add batch correctied assay and batch-corrected RSDQCs
#' exp <- addBatchCorrection(exp)
#'
#' # Show RSDQCs
#' rowData(exp)$rsdqc
#'
#' # Show batch-corrected RSDQCs
#' rowData(exp)$rsdqcCorrected
addBatchCorrection <- function(
        exp, assay = "ratio", qcType = metadata(exp)$QC,
        removeOutliers = TRUE, useWithinBatch = TRUE,
        calculateRSDQC = FALSE
) {
    # Validate if the experiment is correct
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    correctedAssay <- sprintf("%s_corrected", assay)
    if (useWithinBatch) {
        exp <- addWithinBatchCorrectionAssay(
            exp = exp,
            assayName = assay,
            saveAssay = correctedAssay
        )
        assay <- correctedAssay
    }

    exp <- addBatchCorrectionAssay(
        exp = exp,
        assayName = assay,
        saveAssay = correctedAssay,
        removeOutliers = removeOutliers
    )

    if (calculateRSDQC) {
        rowData(exp)$rsdqc <- rsdqc(exp, assay = assay, type = qcType)
        rowData(exp)$rsdqcCorrected <- rsdqc(
            exp = exp, assay = correctedAssay, type = qcType
        )
    }


    # Return the updated experiment
    return(exp)
}
