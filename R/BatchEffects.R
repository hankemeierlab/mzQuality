betweenBatchCorrection <- function(ratio, qcAliquots, qcBatches, aliquots, batches) {

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
        }

        # Multiply the entire batch with the calculated correction factors
        aliqs <- aliquots[batches == batch]
        ratio[, aliqs] <- ratio[, aliqs] * correction_factor
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
addBatchCorrectionAssay <- function(exp, assayName = "ratio", saveAssay = "ratio_corrected", qc = metadata(exp)$QC, removeOutliers = TRUE) {
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
#' @param assay Character string specifying the name of the assay to correct
#'   (default is "ratio").
#' @param saveAssay Character string specifying the name of the new assay to
#'   store corrected values (default is "ratio_corrected").
#' @param qcType Character string specifying the sample type used for QC
#'   (usually "SQC" or "QC"). Defaults to the value stored in
#'   `metadata(exp)$QC`.
#' @param removeOutliers Logical indicating whether to remove extreme QC
#'   outliers at the batch front and back before correction (default is TRUE).
#'
#' @importFrom dplyr group_by row_number pull n mutate arrange filter
#' @export
addWithinBatchCorrectionAssay <- function(exp, assay = "ratio", saveAssay = "ratio_corrected",
                                          qcType = metadata(exp)$QC, removeOutliers = TRUE) {
    x <- lapply(unique(exp$batch), function(j) {
        x <- exp[, exp$batch == j]
        qcs <- x[, x$type == qcType]
        vals <- assay(qcs, assay)


        orderMatrix <- matrix(
            data = rep(x$order[x$type == qcType], nrow(vals)),
            ncol = ncol(qcs),
            byrow = TRUE
        )


        slope <- rowWiseSlope(orderMatrix, vals)
        int <- rowWiseIntercept(orderMatrix, y = vals, slope = slope)

        # Calculate the raw residuals by subtracting the intercept and slope * x-vals
        residuals <- vals - int - slope * orderMatrix

        # Calculate the studentized-residuals (t-distributed error residuals)
        studentResiduals <- studentResiduals(
            residuals = residuals,
            n = rowSums(!is.na(vals))
        )

        outliers <- abs(studentResiduals) > 2
        outliers[!is.finite(outliers)] <- FALSE
        mostlyGood <- colSums(outliers) / nrow(outliers) < 0.7

        if (removeOutliers) {
            # only remove values from the first and last column if they are outliers
            startOutliers <- outliers[, 1]

            N <- ncol(outliers)
            endOutliers <- outliers[, N]

            vals[, 1][startOutliers] <- NA
            vals[, N][endOutliers] <- NA
            vals <- vals[, mostlyGood]

            # Repeat the same for x-values
            orderMatrix[, 1][startOutliers] <- NA
            orderMatrix[, N][endOutliers] <- NA
            orderMatrix <- orderMatrix[, mostlyGood]


            # Recalculate the slope, now with added NA values for outliers
            # These will be ignored due to na.rm = TRUE
            slope <- rowWiseSlope(orderMatrix, vals)
            int <- rowWiseIntercept(orderMatrix, y = vals, slope = slope)
        }

        orderMatrix <- matrix(
            data = rep(seq_len(ncol(x)), nrow(vals)),
            ncol = ncol(x),
            byrow = TRUE
        )
        vals <- orderMatrix * slope + int

        return(
            list(
                Predicted = vals,
                StudentResiduals = studentResiduals,
                Outliers = outliers
            )
        )
    })

    all <- do.call(cbind, lapply(x, `[[`, "Predicted"))
    residualMatrix <- do.call(cbind, lapply(x, `[[`, "StudentResiduals"))
    outlierMatrix <- do.call(cbind, lapply(x, `[[`, "Outliers"))

    dimnames(all) <- dimnames(exp)
    m <- assay(exp, assay) / as.matrix(all)




    m <- do.call(cbind, lapply(unique(exp$batch), function(j) {
        m[, exp$batch == j] * rowMedians(assay(exp[, exp$batch == j & exp$type == "SQC"], assay), na.rm = TRUE)
    }))

    m[!is.finite(m)] <- NA
    m[m < 0] <- NA

    if (removeOutliers) {
        aliqs <- as.data.frame(colData(exp)) %>%
            mutate(aliquot = colnames(exp)) %>%
            arrange(.data$injection_time) %>%
            filter(.data$type == metadata(exp)$QC) %>%
            group_by(.data$batch)

        front <- aliqs %>%
            filter(row_number() == 1) %>%
            pull(.data$aliquot)

        back <- aliqs %>%
            filter(row_number() == n()) %>%
            pull(.data$aliquot)

        m[, front][which(outlierMatrix[, front])] <- NA
        m[, back][which(outlierMatrix[, back])] <- NA
    }


    assay(exp, saveAssay) <- m[, colnames(exp)]
    exp
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
addBatchCorrection <- function(exp, assay = "ratio", qcType = metadata(exp)$QC,
                               removeOutliers = TRUE, useWithinBatch = TRUE) {
    # Validate if the experiment is correct
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    rowData(exp)$rsdqc <- rsdqc(exp, assay = assay, type = qcType)


    correctedAssay <- sprintf("%s_corrected", assay)
    if (useWithinBatch) {
        exp <- addWithinBatchCorrectionAssay(
            exp = exp,
            assay = assay,
            saveAssay = correctedAssay,
            removeOutliers = removeOutliers
        )
        assay <- correctedAssay
    }


    exp <- addBatchCorrectionAssay(exp, assayName = assay, saveAssay = correctedAssay, removeOutliers = removeOutliers)

    rowData(exp)$rsdqcCorrected <- rsdqc(exp, assay = correctedAssay, type = qcType)
    # Return the updated experiment
    return(exp)
}
