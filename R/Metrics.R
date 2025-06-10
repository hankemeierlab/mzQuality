#' @title Calculate the ratio between QCs and SAMPLEs
#' @description Computes the ratio of median values between QC samples and
#'     SAMPLEs for a specified assay in a SummarizedExperiment object.
#' @details This function calculates the ratio of row medians for QC samples
#'     and SAMPLEs in the specified assay. The result is stored in the
#'     `rowData` of the SummarizedExperiment object under the column
#'     `ratioQcSample`.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param assay A string specifying the assay to use for the calculation.
#'     Defaults to `"ratio_corrected"`.
#' @param qcLabel A string specifying the label for QC samples. Defaults to
#'     the `QC` metadata of the experiment.
#' @param sampleLabel A string specifying the label for SAMPLEs. Defaults to
#'     `"SAMPLE"`.
#' @return A SummarizedExperiment object with the calculated ratio added to
#'     the `rowData`.
#' @importFrom SummarizedExperiment rowData<-
#' @export
#' @examples
#' # Example usage:
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#' exp <- ratioQcSample(exp, assay = "ratio")
ratioQcSample <- function(
        exp, assay = "ratio_corrected",
        qcLabel = metadata(exp)$QC, sampleLabel = "SAMPLE"
) {
    # Check if the experiment is valid
    stopifnot(isValidExperiment(exp))

    useExp <- exp[, exp$use]

    if (!all(c(qcLabel, sampleLabel) %in% useExp$type)) {
        return(exp)
    }

    # Subset the experiment with QCs
    qcs <- assay(useExp[, useExp$type == qcLabel], assay)

    # Subset the experiment with SAMPLEs
    samps <- assay(useExp[, useExp$type == sampleLabel], assay)

    # Calculate the Ratio
    ratio <- rowMedians(qcs, na.rm = TRUE) / rowMedians(samps, na.rm = TRUE)
    ratio[!is.finite(ratio)] <- NA
    rowData(exp)$ratioQcSample <- ratio

    # Return the updated experiment
    return(exp)
}

#' @title Identify outliers in QC samples
#' @description Identifies outliers in QC samples for a specified assay in a
#'     SummarizedExperiment object and updates the `colData` with outlier
#'     information.
#' @details This function detects outliers in QC samples based on the
#'     distribution of log-transformed assay values. It uses Rosner's test to
#'     identify outliers and updates the `colData` of the experiment with
#'     boolean columns `use` and `outlier`. Outliers are marked as `FALSE` in
#'     the `use` column and `TRUE` in the `outlier` column.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param assay A string specifying the assay to use for outlier detection.
#'     Defaults to `"ratio"`.
#' @param qcType A string specifying the type of QC samples. Defaults to the
#'     `QC` metadata of the experiment.
#' @return A SummarizedExperiment object with updated `colData` containing
#'     `use` and `outlier` columns.
#' @importFrom matrixStats colMedians
#' @importFrom EnvStats rosnerTest
#' @importFrom SummarizedExperiment colData assay
#' @export
#' @examples
#' # Example usage:
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#' exp <- identifyOutliers(exp, assay = "ratio")
identifyOutliers <- function(exp, assay = "ratio", qcType = metadata(exp)$QC) {
    stopifnot(isValidExperiment(exp))

    if ("use" %in% colnames(colData(exp))) {
        exp$use <- colData(exp)$use
    } else {
        exp$use <- TRUE
    }

    if ("outlier" %in% colnames(colData(exp))) {
        exp$outlier <- colData(exp)$outlier
    } else {
        exp$outlier <- FALSE
    }

    qc_se <- exp[, exp$type == qcType]

    if (nrow(qc_se) <= 3) {
        return(exp)
    }

    medians <- colMedians(log10(assay(qc_se, assay)), na.rm = TRUE)
    medians <- medians[is.finite(medians)]
    non_na <- sum(!is.na(medians))

    if (non_na < 3) {
        return(exp)
    }

    batches <- length(unique(exp$batch))
    k <- ifelse(non_na > batches, batches, non_na - 2)
    test <- EnvStats::rosnerTest(medians, k = k, alpha = 0.05, warn = FALSE)
    outliers <- test$all.stats$Obs.Num[test$all.stats$Outlier]

    idx <- match(names(medians[outliers]), colnames(exp))
    if (length(idx) > 0) {
        exp$use[idx] <- FALSE
        exp$outlier[idx] <- TRUE
    }

    return(exp)
}

#' @title Identify misinjections in Samples
#' @description Identifies misinjections in SAMPLEs for a specified assay in
#'     a SummarizedExperiment object.
#' @details This function is a wrapper around `identifyOutliers()` to detect
#'     misinjections in SAMPLEs based on the distribution of assay values.
#'     It updates the `colData` of the experiment with boolean columns `use`
#'     and `outlier`.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param assay A string specifying the assay to use for outlier detection.
#'     Defaults to the `secondary` metadata of the experiment.
#' @param type A string specifying the type of samples to analyze. Defaults
#'     to `"SAMPLE"`.
#' @return A SummarizedExperiment object with updated `colData` containing
#'     `use` and `outlier` columns.
#' @examples
#' # Example usage:
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#' exp <- identifyMisInjections(exp, assay = "area_is")
#' @export
identifyMisInjections <- function(
        exp, assay = metadata(exp)$secondary, type = "SAMPLE"
) {
    identifyOutliers(
        exp = exp,
        assay = metadata(exp)$secondary,
        qcType = type
    )
}

#' @title Calculate the presence of a specific type in the experiment
#' @description Computes the proportion of non-missing values for a specific
#'     type (e.g., QC) in the assay data of a SummarizedExperiment object.
#' @details This function calculates the proportion of non-missing values for
#'     a specified type (e.g., QC) in the assay data. The result is stored in
#'     the `rowData` of the experiment under a column named
#'     `<type>Presence`.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param area A string specifying the assay to use for the calculation.
#'     Defaults to the `primary` metadata of the experiment.
#' @param type A string specifying the type of samples to analyze. Defaults
#'     to the `QC` metadata of the experiment.
#' @return A SummarizedExperiment object with the calculated presence added
#'     to the `rowData`.
typePresence <- function(
        exp, area = metadata(exp)$primary,
        type = metadata(exp)$QC
) {
    useExp <- exp[, exp$use]
    areas <- assay(useExp, area)
    cols <- useExp$type == type
    rowData(exp)[[sprintf("%sPresence", type)]] <-
        rowSums(!is.na(areas[, cols, drop = FALSE])) / sum(cols)
    return(exp)
}

#' @title Calculate the median sample area for each compound
#' @description Computes the median area for each compound in SAMPLEs for a
#'     specified assay in a SummarizedExperiment object.
#' @details This function calculates the median area for each compound in
#'     SAMPLEs based on the specified assay. The result is stored in the
#'     `rowData` of the experiment under the column `MedianSampleArea`.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param type A string specifying the type of samples to analyze. Defaults
#'     to `"SAMPLE"`.
#' @param assayName A string specifying the assay to use for the calculation.
#'     Defaults to the `primary` metadata of the experiment.
#' @return A SummarizedExperiment object with the calculated median sample
#'     area added to the `rowData`.
#' @importFrom stats median
medianSampleArea <- function(
        exp, type = "SAMPLE", assayName = metadata(exp)$primary
) {
    m <- assay(exp[, exp$type == type], assayName)
    medians <- unlist(apply(m, 1, median, na.rm = TRUE))

    rowData(exp)$MedianSampleArea <- medians

    return(exp)
}

#' @title Calculate compound usability metrics
#' @description Evaluates the usability of compounds in a SummarizedExperiment
#'     object based on RSDQC, background signal, and QC presence thresholds.
#' @details This function calculates whether each compound meets the
#'     predefined thresholds for RSDQC, background signal, and QC presence.
#'     The result is stored in the `rowData` of the experiment under the
#'     column `use`.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param backgroundPercentage Numeric, the maximum allowable percentage of
#'     background signal. Defaults to 40.
#' @param qcPercentage Numeric, the minimum percentage of QC presence
#'     required. Defaults to 80.
#' @param nonReportableRSD Numeric, the maximum allowable RSDQC for compounds
#'     to be considered usable. Defaults to 30.
#' @return A SummarizedExperiment object with the usability metrics added to
#'     the `rowData`.
#' @noRd
.calculateMetrics <- function(
        exp, backgroundPercentage = 40,
        qcPercentage = 80, nonReportableRSD = 30
) {
    df <- rowData(exp)

    rowData(exp)$rsdqc <- rsdqc(
        exp = exp,
        assay = "ratio",
        type = metadata(exp)$QC
    )

    rowData(exp)$rsdqcCorrected <- rsdqc(
        exp = exp,
        assay = "ratio_corrected",
        type = metadata(exp)$QC
    )

    rsdThreshold <- df$rsdqcCorrected <= nonReportableRSD
    column <- sprintf("%sPresence", metadata(exp)$QC)
    qcPresenceThreshold <- df[, column] * 100 >= qcPercentage

    lowBackground <- df$backgroundSignal * 100 <= backgroundPercentage
    noBackground <- is.na(df$backgroundSignal)
    backgroundThreshold <- lowBackground | noBackground

    df$use <- rsdThreshold & backgroundThreshold & qcPresenceThreshold
    df$use[is.na(df$use)] <- FALSE

    rowData(exp) <- df
    return(exp)
}

#' @title Perform the main analysis of mzQuality
#' @description Executes the main analysis pipeline for mzQuality, including
#'     batch correction, RSDQC calculation, background signal estimation, and
#'     internal standard suggestion.
#' @details This function validates the input SummarizedExperiment object and
#'     performs a series of analyses, including ratio calculation, batch
#'     correction, background signal estimation, and RSDQC calculation. If
#'     concentration data is available, it also calculates concentrations and
#'     applies additional corrections. The function can optionally filter out
#'     bad compounds and outliers.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param removeBadCompounds Logical, whether to remove compounds that do not
#'     meet usability thresholds. Defaults to `FALSE`.
#' @param removeOutliers Logical, whether to remove outliers during batch
#'     correction. Defaults to `FALSE`.
#' @param useWithinBatch Logical, whether to use within-batch correction.
#'     Defaults to `FALSE`.
#' @param effectNaAsZero Logical, whether to treat missing values as zero
#'     during background signal estimation. Defaults to `FALSE`.
#' @param concentrationType A string specifying the concentration type to use
#'     for calculations. Defaults to `"ACAL"`.
#' @param backgroundPercentage Numeric, the maximum allowable percentage of
#'     background signal. Defaults to 40.
#' @param qcPercentage Numeric, the minimum percentage of QC presence
#'     required. Defaults to 80.
#' @param nonReportableRSD Numeric, the maximum allowable RSDQC for compounds
#'     to be considered usable. Defaults to 30.
#' @return A SummarizedExperiment object with added assays and metadata
#'     reflecting the analysis results.
#' @examples
#' # Read data
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#'
#' # Do Analysis
#' exp <- doAnalysis(exp)
#' @export
doAnalysis <- function(
        exp, removeBadCompounds = FALSE, removeOutliers = FALSE,
        useWithinBatch = FALSE, effectNaAsZero = FALSE,
        concentrationType = "ACAL", backgroundPercentage = 40,
        qcPercentage = 80, nonReportableRSD = 30
) {
    stopifnot(isValidExperiment(exp))

    exp <- calculateRatio(exp) %>%
        addBatchCorrection(
            removeOutliers = removeOutliers,
            useWithinBatch = useWithinBatch,
            calculateRSDQC = TRUE
        ) %>%
        backgroundSignals(NaAsZero = effectNaAsZero) %>%
        ratioQcSample() %>%
        typePresence() %>%
        medianSampleArea() %>%
        getSuggestedInternalStandards(
            removeOutliers = removeOutliers,
            useWithinBatch = useWithinBatch
        )

    if ("concentration" %in% assayNames(exp)) {
        assayName <- sprintf("%s_concentration", metadata(exp)$concentration)
        exp <- exp %>%
            calculateConcentrations() %>%
            addBatchCorrection(assay = assayName)
    }

    exp <- .calculateMetrics(
        exp = exp, backgroundPercentage = backgroundPercentage,
        qcPercentage = qcPercentage, nonReportableRSD = nonReportableRSD
    )

    if (removeBadCompounds) {
        exp <- exp[rowData(exp)$use, ]
    }
    return(exp)
}

#' @title Suggest internal standards for compounds
#' @description Identifies the internal standard that minimizes the
#'     batch-corrected RSDQC for each compound in a SummarizedExperiment
#'     object.
#' @details This function evaluates all available internal standards and
#'     suggests the one that results in the lowest batch-corrected RSDQC for
#'     each compound. The results are stored in the `rowData` of the
#'     experiment under the columns `suggestedIS` and `suggestedRSDQC`.
#' @param exp A SummarizedExperiment object containing the experimental
#'     data.
#' @param secondaryAssay A string specifying the assay to use for the
#'    calculation. Defaults to the `secondary` metadata of the experiment.
#' @param removeOutliers Logical, whether to remove outliers during batch
#'     correction. Defaults to `TRUE`.
#' @param useWithinBatch Logical, whether to use within-batch correction.
#'     Defaults to `TRUE`.
#' @return A SummarizedExperiment object with the suggested internal standards
#'     and their corresponding RSDQC values added to the `rowData`.
#' @export
#' @examples
#' # Example usage:
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#' exp <- doAnalysis(exp, removeOutliers = TRUE)
#' exp <- getSuggestedInternalStandards(exp, removeOutliers = TRUE)
getSuggestedInternalStandards <- function(
        exp, secondaryAssay = metadata(exp)$secondary,
        removeOutliers = TRUE, useWithinBatch = TRUE
) {
    if (!metadata(exp)$hasIS) {
        return(exp)
    }

    rowData(exp)$compound <- rownames(exp)
    comp_is <- unique(rowData(exp)$compound_is)
    df <- expand.grid(
        compound_is = comp_is,
        compound = rownames(exp)
    )
    comp_row <- match(df$compound, rownames(exp))
    testExp <- exp[comp_row, ]
    rownames(testExp) <- seq_len(nrow(testExp))

    comp_is_row <- match(df$compound_is, rowData(testExp)$compound_is)

    m <- assay(testExp, secondaryAssay)[comp_is_row, , drop = FALSE]
    assay(testExp, secondaryAssay, withDimnames = FALSE) <- m

    rowData(testExp)$compound_is <- as.character(df$compound_is)

    metrics <- testExp %>%
        calculateRatio() %>%
        addBatchCorrection(
            removeOutliers = removeOutliers,
            useWithinBatch = useWithinBatch,
            calculateRSDQC = TRUE
        )

    exp <- .addSuggestedMetrics(exp, rowData(metrics))

    return(exp)
}

#' @title Add suggested metrics to the experiment
#' @description Adds suggested internal standards and their corresponding
#' RSDQC values to the `rowData` of the provided SummarizedExperiment
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param metrics A data frame containing the suggested internal standards and
#' their corresponding RSDQC values.
#' @importFrom dplyr group_by summarise mutate %>%
#' @importFrom SummarizedExperiment rowData<-
#' @noRd
.addSuggestedMetrics <- function(exp, metrics){
    if (nrow(metrics) == 0) {
        return(exp)
    }

    metrics <- as.data.frame(metrics) %>%
        group_by(.data$compound) %>%
        mutate(hasRSDQC = any(!is.na(.data$rsdqcCorrected))) %>%
        summarise(
            suggestedRSDQC = ifelse(
                .data$hasRSDQC[1],
                min(.data$rsdqcCorrected, na.rm = TRUE),
                NA
            ),
            suggestedIS = ifelse(
                .data$hasRSDQC[1],
                .data$compound_is[which.min(.data$rsdqcCorrected)],
                NA
            )
        )

    rowData(exp)$suggestedIS <- metrics$suggestedIS
    rowData(exp)$suggestedRSDQC <- metrics$suggestedRSDQC
    return(exp)
}
