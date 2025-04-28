#' @title Calculate the ratio between QCs and SAMPLEs
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp
#' @param assay
#' @param qcLabel
#' @param sampleLabel
#' @importFrom SummarizedExperiment rowData<-
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Calculate the ratio between QCs and Samples of the batch-corrected ratio
#' ratioQcSample(exp)
ratioQcSample <- function(exp, assay = "ratio_corrected",
                          qcLabel = metadata(exp)$QC, sampleLabel = "SAMPLE") {
    # Check if the experiment is valid
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

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

#' @title Identifier outliers in QC
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment Object
#' @param assay Which assay should be used for outlier detection
#' @importFrom matrixStats colMedians
#' @returns SummarizedExperiment with added column called "Use" with boolean
#' values suggesting if it should be used in analysis
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Identify aliquot outliers and add to colData
#' exp <- identifyOutliers(exp, assay = "Ratio")
identifyOutliers <- function(exp, assay = "ratio", qcType = metadata(exp)$QC) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    if ("use" %in% colnames(colData(exp))) {
        exp$use <- colData(exp)$use
        exp$outlier <- colData(exp)$outlier
    } else {
        exp$use <- TRUE
        exp$outlier <- FALSE
    }


    if (!assay %in% assayNames(exp)) {
        return(exp)
    }
    if (!qcType %in% exp$type) {
        return(exp)
    }


    qc_se <- exp[, exp$type == qcType]

    if (nrow(qc_se) <= 3) {
        return(exp)
    }

    medians <- colMedians(log10(assay(qc_se, assay)), na.rm = TRUE)

    idx <- is.finite(medians)

    medians <- medians[idx]
    non_na <- sum(!is.na(medians))


    if (non_na < 3) {
        return(exp)
    }


    batches <- length(unique(exp$batch))

    k <- ifelse(non_na > batches, batches, non_na - 2)
    test <- EnvStats::rosnerTest(medians, k = k, alpha = 0.05, warn = FALSE)
    outliers <- test$all.stats$Obs.Num[test$all.stats$Outlier]

    idx <- match(names(medians[outliers]), colnames(exp))
    exp$use[idx] <- FALSE
    exp$outlier[idx] <- TRUE
    exp
}

#' @title Identifier misinjections in Samples
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment Object
#' @param assay Which assay should be used for outlier detection
#' @param type
#' @importFrom matrixStats colMedians
#' @importFrom EnvStats rosnerTest
#' @importFrom SummarizedExperiment colData assay
#' @returns SummarizedExperiment with added column called "Use" with boolean
#' values suggesting if it should be used in analysis
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Identify aliquot outliers and add to colData
#' exp <- identifyOutliers(exp, assay = "Ratio")
identifyMisInjections <- function(exp, assay = metadata(exp)$secondary, type = "SAMPLE") {

    identifyOutliers(
        exp = exp,
        assay = metadata(exp)$secondary,
        qcType = type
    )
}

typePresence <- function(exp, area = metadata(exp)$primary,
                         type = metadata(exp)$QC) {
    useExp <- exp[, exp$use]
    areas <- assay(useExp, area)
    cols <- useExp$type == type
    rowData(exp)[[sprintf("%sPresence", type)]] <- rowSums(!is.na(areas[, cols, drop = FALSE])) / sum(cols)
    return(exp)
}

medianSampleArea <- function(exp, type = "SAMPLE", assayName = metadata(exp)$primary) {
    m <- assay(exp[, exp$type == type], assayName)
    medians <- unlist(apply(m, 1, median, na.rm = TRUE))

    rowData(exp)$MedianSampleArea <- medians

    return(exp)
}

#' @title Perform the main analysis of mzQuality
#' @description This function is the main analysis function of mzQuality. It
#' starts by validating the SummarizedExperiment object is valid. Next, it will
#' calculate the batch-corrected RSDQC, background signals, retention time
#' shifts and lastly calculates a suggested internal standard that results in
#' the lowest batch-corrected RSDQC.
#' @details placeholder
#' @returns SummarizedExperiment with added slots and assays with
#' batch-corrected values, including RSDQC, background signal, rt shift and
#' suggested internal standards.
#' @param exp SummarizedExperiment Object
#' @param aliquots
#' @param doAll
#' @importFrom dplyr %>%
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Perform analysis
#' exp <- doAnalysis(exp)
#' exp
doAnalysis <- function(exp,
                       aliquots = colnames(exp),
                       doAll = FALSE,
                       removeBadCompounds = FALSE,
                       removeOutliers = FALSE,
                       useWithinBatch = FALSE,
                       effectNaAsZero = FALSE,
                       concentrationType = "ACAL",
                       backgroundPercentage = 40,
                       qcPercentage = 80,
                       nonReportableRSD = 30) {
    message("Start Analysis")

    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    doAliquots <- length(aliquots) != ncol(exp)

    exp <- exp[, aliquots]

    if (doAliquots | doAll | !"rsdqc" %in% colnames(rowData(exp))) {
        exp <- calculateRatio(exp) %>%
            addBatchCorrection(removeOutliers = removeOutliers, useWithinBatch = useWithinBatch) %>%
            backgroundSignals(NaAsZero = effectNaAsZero) %>%
            matrixEffect() %>%
            ratioQcSample() %>%
            typePresence() %>%
            medianSampleArea() %>%
            suggestedInternalStandards(
                removeOutliers = removeOutliers,
                useWithinBatch = useWithinBatch
            )
    }

    if ("concentration" %in% assayNames(exp)) {
        exp <- exp %>%
            calculateConcentrations(type = concentrationType) %>%
            addBatchCorrectionAssay(assay = "concentration") %>%
            carryOverEffect() %>%
            blankLimits()
    }

    rsdThreshold <- rowData(exp)$rsdqcCorrected <= nonReportableRSD

    qcPresenceThreshold <- rowData(exp)[, sprintf("%sPresence", metadata(exp)$QC)] * 100 >= qcPercentage
    backgroundThreshold <- rowData(exp)$backgroundSignal * 100 <= backgroundPercentage | !is.finite(rowData(exp)$backgroundSignal)
    rowData(exp)$qcPresenceThreshold <- qcPresenceThreshold

    rowData(exp)$use <- rsdThreshold & backgroundThreshold & qcPresenceThreshold
    rowData(exp)$use[is.na(rowData(exp)$use)] <- FALSE

    if (removeBadCompounds) {
        exp <- exp[rowData(exp)$use, ]
    }
    message("Finished Analysis")



    exp
}

suggestedInternalStandards <- function(experiment, removeOutliers, useWithinBatch) {
    if (!metadata(experiment)$hasIS) {
        return(experiment)
    }

    rowData(experiment)$compound <- rownames(experiment)
    comp_is <- unique(rowData(experiment)$compound_is)
    df <- expand.grid(
        compound_is = comp_is,
        compound = rownames(experiment)
    )
    comp_row <- match(df$compound, rownames(experiment))
    exp <- experiment[comp_row, ]
    rownames(exp) <- seq_len(nrow(exp))

    comp_is_row <- match(df$compound_is, rowData(exp)$compound_is)

    assay(exp, "area_is", withDimnames = FALSE) <- assay(exp, "area_is")[comp_is_row, , drop = FALSE]
    rowData(exp)$compound_is <- df$compound_is

    testExp <- calculateRatio(exp) %>%
        addBatchCorrection(removeOutliers = removeOutliers, useWithinBatch = useWithinBatch)

    x <- as.data.frame(rowData(testExp))

    rowData(experiment)$suggestedIS <- unlist(lapply(split(x, x$compound), function(z) {
        if (all(is.na(z$rsdqcCorrected))) {
            return(NA)
        }
        as.character(z$compound_is[which.min(z$rsdqcCorrected)])
    }))

    rowData(experiment)$suggestedRSDQC <- unlist(lapply(split(x, x$compound), function(z) {
        if (all(is.na(z$rsdqcCorrected))) {
            return(NA)
        }
        min(z$rsdqcCorrected, na.rm = TRUE)
    }))
    experiment
}
