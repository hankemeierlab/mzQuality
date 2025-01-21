#' @title Calculate Carry-over Effect
#' @description The carry-over effect represents a percentage of the compound
#' concentration that has been measured in a blank sample following another
#' sample. Ideally this value should be as low as possible. In order to run
#' this function, the calculations for all samples must be calculated using
#' `calculateConcentrations()`.
#' @param exp A SummarizedExperiment object
#' @param type Name of the sample type that represents a blank sample,
#' defaults to 'BLANK'.
#' @returns SummarizedExperiment with an added column in the rowData slot
#' called "CarryOver"
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
#'
#' # Add concentrations
#' concentrations <- read.delim(system.file(
#'     package = "mzQuality",
#'     "concentrations.txt"
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
carryOverEffect <- function(exp, type = "PROC", assay = "concentration",
                            method = c("previous", "calibration")[1]) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    rowData(exp)$carryOver <- NA

    blanks <- which(exp$type == type)

    if (method == "previous") {
      other <- which(exp$type != type)
      blanks <- blanks[which((blanks - 1) %in% other)]

      if (length(blanks) > 0 & length(other) > 0 & assay %in% assayNames(exp)) {
        df <- do.call(cbind, lapply(blanks, function(blank) {
          assay1 <- assay(exp[, blank], assay)
          assay2 <- assay(exp[, blank - 1], assay)
          return(assay1 / assay2)
        }))
        from <- colnames(exp)[blanks - 1]
        to <- colnames(exp)[blanks]
        colnames(df) <- glue::glue("{from}_{to}", from = from, to = to)
        rowData(exp)$carryOver <- df
      }
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
#' @returns SummarizedExperiment with 3 added columns in the rowData slot
#' called "LoB" (Limit of Blank), "LoD" (Limit of Detection), and "LoQ"
#' (Limit of Quantification)
#' @export
#' @examples
#' # Read example dataset
#' data <- read.delim(system.file(package = "mzQuality", "dataset.txt"))
#'
#' # Construct experiment
#' exp <- buildExperiment(
#'     data,
#'     rowIndex = "Compound",
#'     colIndex = "Aliquot",
#'     primaryAssay = "Area",
#'     secondaryAssay = "Area_is"
#' )
#'
#' # Add concentrations
#' concentrations <- read.delim(system.file(
#'     package = "mzQuality",
#'     "concentrations.txt"
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
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }
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
#' @returns SummarizedExperiment with an added "backgroundSignal" column in
#' the rowData slot
#' @export
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
#'
#' # Perform the calculation
#' backgroundSignals(exp)
backgroundSignals <- function(exp, type = "BLANK", NaAsZero = FALSE) {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    rowData(exp)$backgroundSignal <- NA
    if (type %in% toupper(exp$type)) {
      rowData(exp)$backgroundSignal <- calculateEffect(exp = exp, type = type, NaAsZero = NaAsZero)
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
#' @examples
#' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
#'
#' # Calculate the Blank effect
#' calculateEffect(exp, type = "BLANK")
#'
#' # Calculate the PROC effect
#' calculateEffect(exp, type = "PROC")
calculateEffect <- function(exp, assay = metadata(exp)$primary,
                            type = "BLANK", sampleLabel = "SAMPLE", NaAsZero = FALSE) {

    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

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

#' @title Calculate the matrix effect and ion suppression
#' @param exp SummarizedExperiment
#' @importFrom matrixStats rowMedians
#' @export
#' @examples
#' #' # Read the example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))
#'
#' matrixEffect(exp)
matrixEffect <- function(exp, is_assay = metadata(exp)$secondary,
                         sampleLabel = "SAMPLE", procBlankLabel = "PROC"){

  rowData(exp)$matrixEffectFactor <- NA
  if (metadata(exp)$hasIS) {
      useExp <- exp[, exp$use]
      if (all(c(sampleLabel, procBlankLabel) %in% useExp$type)) {

          samps <- rowMeans(assay(useExp[, useExp$type == sampleLabel], is_assay), na.rm = TRUE)
          procs <- rowMeans(assay(useExp[, useExp$type == procBlankLabel], is_assay), na.rm = TRUE)

          # rowData(exp)$matrixEffect <- round(rowMedians(samps) / rowMedians(procs), 3)
          # rowData(exp)$ionSuppresion <- round(100 - rowData(exp)$matrixEffect, 3)

          rowData(exp)$matrixEffectFactor <- round((procs - samps) / procs, 3)
      }
  }

  return(exp)
}
