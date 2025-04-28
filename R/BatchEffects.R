#' @title Calculate the batch effection corrections
#' @param exp SummarizedExperiment Object
#' @importFrom matrixStats rowMedians
#' @returns SummarizedExperiment with batch-effect corrected ratios in the
#' assay slot with name "ratio Corrected.
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


    if (qc %in% exp$type & "batch" %in% colnames(colData(exp))) {
        # subset the experiment with only QCs
        qcExp <- exp[, exp$type == qc]

        # Calculate the median ratio per compound, across all batches
        compound_qc_ratio_median <- rowMedians(ratios[, colnames(qcExp), drop = FALSE], na.rm = TRUE)

        # For each batch, calculate the correction factor and multiply the
        # ratio to obtain the ratio Corrected factors
        df <- do.call(cbind, lapply(unique(exp$batch), function(batch) {
            correction_factor <- 1

            # Check if the batch is present in the QC-subsetted experiment
            if (batch %in% qcExp$batch) {
                # Subset the QC experiment for the current batch
                qcBatchExp <- qcExp[, qcExp$batch == batch]

                # Calculate the median ratio for this batch
                med_ratio_qc <- rowMedians(ratios[, colnames(qcBatchExp), drop = FALSE], na.rm = TRUE)

                # Calculate the QC correction factor by dividing the overall median
                # by this batch median
                correction_factor <- compound_qc_ratio_median / med_ratio_qc
            }
            # Multiply the entire batch with the calculated correction factors
            assay(exp[, exp$batch == batch], assayName) * correction_factor
        }))

        # Ensure that the right column names are used, in the right order
        df <- df[, colnames(exp), drop = FALSE]

        # Ensure that the dimensions of the data.frame are equal to the experiments
        dimnames(df) <- dimnames(exp)

        # Set the ratio corrected values in the assay
        assay(exp, saveAssay) <- df
    }


    m <- assay(exp, saveAssay)
    m[!is.finite(m)] <- NA
    m[, !exp$use] <- NA
    assay(exp, saveAssay) <- m

    # Return the experiment
    return(exp)
}



#' @title Within-Batch correction
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param saveAssay
#' @param qcType
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
