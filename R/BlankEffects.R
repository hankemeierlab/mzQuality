#' @title Calculate Carry-over Effect
#' @description The carry-over effect represents a percentage of the compound
#' concentration that has been measured in a blank sample following another
#' sample. Ideally this value should be as low as possible. In order to run
#' this function, the calculations for all samples must be calculated using
#' `calculateConcentrations()`.
#' @param exp A SummarizedExperiment object
#' @param type Name of the sample type that represents a blank sample,
#' defaults to 'BLANK'.
#' @param assay Name of the assay to use for the calculation, defaults to
#' 'concentration'.
#' @returns SummarizedExperiment with an added column in the rowData slot
#' called "CarryOver"
#' @export
#' @examples
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#'
#' # Add concentrations
#' concentrations <- read.delim(system.file(
#'     package = "mzQuality",
#'     "extdata/concentrations.tsv"
#' ))
#' exp <- addConcentrations(exp, concentrations)
#'
#' # Perform analysis
#' exp <- doAnalysis(exp)
#'
#' # Model concentrations
#' exp <- calculateConcentrations(exp)
#'
#' # Calculate carryOverEffect
#' exp <- carryOverEffect(exp)
#' exp
carryOverEffect <- function(exp, type = "PROC", assay = "concentration") {
    stopifnot(isValidExperiment(exp))

    rowData(exp)$carryOver <- NA

    blanks <- which(exp$type == type)

    other <- which(exp$type != type)
    blanks <- blanks[which((blanks - 1) %in% other)]
    canCompare <- length(blanks) > 0 & length(other) > 0

    if (canCompare & assay %in% assayNames(exp)) {
        df <- do.call(cbind, lapply(blanks, function(blank) {
            assay1 <- assay(exp[, blank], assay)
            assay2 <- assay(exp[, blank - 1], assay)
            return(assay1 / assay2)
        }))
        from <- colnames(exp)[blanks - 1]
        to <- colnames(exp)[blanks]

        colnames(df) <- sprintf("%s_%s", from, to)
        rowData(exp)$carryOver <- df
    }
    return(exp)
}

#' @title Add limits of Blank, Detection and Quantification with Concentrations
#' @description The Limit of Blank (LoB), -Detection (LoD) and -Quantification
#' (LoQ) are metrics used to determine if the calculated absolute
#' quantification is reliable. Compounds with calculated concentrations that
#' approach these limits might be unreliable. The LoB represent the 95%
#' percentile of the Blank samples, the LoD represents the limit where a
#' compound can reliably be detected. The LoQ represents the limit where a
#' compound can reliably be quantified. The results of the LoB, LoD, and LoQ
#' are stored in the rowData of the experiment.
#' @param exp A SummarizedExperiment object
#' @param type Name of the sample type representing processed blanks,
#' defaults to 'PROC'
#' @param assay Name of the assay to use for the calculation, defaults to
#' 'concentration'
#' @returns SummarizedExperiment with 3 added columns in the rowData slot
#' called "LoB" (Limit of Blank), "LoD" (Limit of Detection), and "LoQ"
#' (Limit of Quantification)
#' @export
#' @examples
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#'
#' # Add concentrations
#' concentrations <- read.delim(system.file(
#'     package = "mzQuality",
#'     "extdata/concentrations.tsv"
#' ))
#' exp <- addConcentrations(exp, concentrations)
#'
#' # Perform analysis
#' exp <- doAnalysis(exp)
#'
#' # Model concentrations
#' exp <- calculateConcentrations(exp)
#'
#' # Calculate limits of blank samples
#' exp <- blankLimits(exp)
#' exp
blankLimits <- function(exp, type = "PROC", assay = "concentration") {
    stopifnot(isValidExperiment(exp))
    rowData(exp)$lob <- NA
    rowData(exp)$lod <- NA
    rowData(exp)$loq <- NA

    if (assay %in% assayNames(exp) & type %in% exp$type) {
        m <- rowMeans(assay(exp[, exp$type == type], assay), na.rm = TRUE)
        s <- rowSds(assay(exp[, exp$type == type], assay), na.rm = TRUE)

        rowData(exp)$lob <- round(m + 1.645 * s, 3)
        rowData(exp)$lod <- round(m + 3.3 * s, 3)
        rowData(exp)$loq <- round(m + 10 * s, 3)
    }

    return(exp)
}

#' @title Calculate background signals in compounds
#' @description Background signal (also called blank effect) is the relative
#' intensity measured in blank samples compared to study samples. Here, the
#' background signal is calculated for the sample with the given blank type.
#' This is done by dividing the means blank signal by the median sample signal.
#' If this value exceeds 0.4, it is generally cosidered having a poor
#' signal-to-noise (SNR) ratio. Consider removing compounds with a poor SNR.
#' @param exp SummarizedExperiment Object
#' @param type Name of the sample type representing blank samples.
#' @param NaAsZero Logical, if TRUE, all NA values are replaced with 0 before
#' calculating the background signal.
#' @returns SummarizedExperiment with an added "backgroundSignal" column in
#' the rowData slot
#' @export
#' @examples
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#'
#' # Perform the calculation
#' backgroundSignals(exp)
backgroundSignals <- function(exp, type = "BLANK", NaAsZero = FALSE) {
    stopifnot(isValidExperiment(exp))

    rowData(exp)$backgroundSignal <- NA
    if (type %in% toupper(exp$type)) {
        rowData(exp)$backgroundSignal <- .calculateEffect(
            exp = exp, type = type, NaAsZero = NaAsZero
        )
    }


    return(exp)
}

#' @title Helper function for calculating effects compared to the sample
#' @param exp SummarizedExperiment Object
#' @param type Sample type to test against a sample
#' @importFrom matrixStats rowMedians
#' @returns A vector of blank effects with equal length to the compounds in
#' the provided SummarizedExperiment
#' @noRd
.calculateEffect <- function(
        exp, assay = metadata(exp)$primary,
        type = "BLANK", sampleLabel = "SAMPLE",
        NaAsZero = FALSE
) {
    stopifnot(isValidExperiment(exp))

    exp <- exp[, exp$use]

    effect <- NA
    if (all(c(sampleLabel, type) %in% exp$type)) {
        blanks <- assay(exp[, exp$type == type], assay)

        if (NaAsZero) {
            blanks[is.na(blanks)] <- 0
        }

        means <- rowMeans(blanks, na.rm = TRUE)

        samples <- assay(exp[, exp$type == sampleLabel], assay)

        if (NaAsZero) {
            samples[is.na(samples)] <- 0
        }
        medians <- rowMedians(samples, na.rm = TRUE)

        effect <- round(means / medians, 5)
    }
    effect[!is.finite(effect)] <- NA
    return(effect)
}

#' @title Calculate Matrix Effect and Ion Suppression in LC-MS Analysis
#' @description Quantifies the matrix effect factor by comparing internal
#'   standard signal intensities between study samples and processed blank
#'   samples.
#' @details Matrix effects in LC-MS analysis occur when components of the sample
#'   matrix enhance or suppress the ionization of target analytes, affecting
#'   quantification accuracy. This function calculates a matrix effect factor
#'   using the formula: (proc_blank_IS - sample_IS) / proc_blank_IS, where IS
#'   refers to internal standard signals.
#'
#'   A matrix effect factor near zero indicates minimal matrix effect, while
#'   negative values indicate ion enhancement and positive values indicate ion
#'   suppression. Values are stored in the rowData of the experiment.
#'
#'   The function only works if internal standards are present in the dataset
#'   (metadata(exp)$hasIS must be TRUE).
#'
#' @param exp A SummarizedExperiment object containing metabolomics data
#' @param is_assay The assay containing internal standard data, defaults to
#'   metadata(exp)$secondary
#' @param sampleLabel Character string indicating the sample type identifier for
#'   study samples, defaults to "SAMPLE"
#' @param procBlankLabel Character string indicating the sample type identifier
#'   for processed blank samples, defaults to "PROC"
#' @return A SummarizedExperiment object with an added "matrixEffectFactor"
#'   column in the rowData slot, containing values between -1 and 1
#' @importFrom matrixStats rowMedians
#' @export
#' @examples
#' path <- system.file("extdata", "example.tsv", package = "mzQuality")
#' exp <- buildExperiment(readData(path))
#'
#' # Calculate matrix effect factors
#' exp <- matrixEffect(exp)
matrixEffect <- function(
        exp, is_assay = metadata(exp)$secondary,
        sampleLabel = "SAMPLE", procBlankLabel = "PROC"
) {

    rowData(exp)$matrixEffectFactor <- NA
    if (metadata(exp)$hasIS) {
        useExp <- exp[, exp$use]
        if (all(c(sampleLabel, procBlankLabel) %in% useExp$type)) {
            mat <- assay(useExp[, useExp$type == sampleLabel], is_assay)
            samps <- rowMeans(mat, na.rm = TRUE)

            mat <- assay(useExp[, useExp$type == procBlankLabel], is_assay)
            procs <- rowMeans(mat, na.rm = TRUE)

            rowData(exp)$matrixEffectFactor <- round((procs - samps) / procs, 3)
        }
    }

    return(exp)
}
