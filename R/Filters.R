#' @title Filter compounds that are equal to given internal standards
#' @description In a scenario where internal standards were measured as
#' compounds, this function aids in removing these from the
#' SummarizedExperiment before doing analysis. It filter by both the `tag`
#' argument and by checking if the A/IS ratio equals 1.
#' @details placeholder
#' @returns SummarizedExperiment with removed compounds based on the filter
#' provided
#' @param exp A SummarizedExperiment object
#' @param tag Tag to search for in compound names that indicates an internal
#' standard. Defaults to 'ISTD'
#' Regular expressions are supported.
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
#' # Filter internal standards from compounds
#' filterISTD(exp, tag = "ISTD")
#'
#' # Filter all standards from compounds
#' filterSTD(exp, tag = "STD")
filterISTD <- function(exp, tag = "STD") {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }

    exp <- exp[grep(tag, rownames(exp), invert = TRUE, ignore.case = TRUE), ]
    if (metadata(exp)$hasIS) {
        a <- assay(exp, metadata(exp)$primary)
        b <- assay(exp, metadata(exp)$secondary)
        exp <- exp[rowSums(a / b == 1 | is.na(a / b), na.rm = TRUE) != ncol(exp), ]
    }
    exp
}

#' @title Filter SST Aliquots from the SummarizedExperiment
#' @description System Sustainability Tests (SST) samples are often executed by
#' the technician for control during measurements. However, they are often
#' unwanted in analysis of measurements. This function helps to remove these
#' SSTs from the SummarizedExperiment object
#' @details placeholder
#' @returns SummarizedExperiment with removed aliquots based on the filter
#' provided
#' @param exp A SummarizedExperiment object
#' @param tag Tag to search for in aliquot types that indicate a SST.
#' Defaults to 'SST'
#' Regular expressions are supported.
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
#' # Filter SST aliquots
#' filterSST(exp, tag = "SST")
filterSST <- function(exp, tag = "SST") {
    if (!validateExperiment(exp)) {
        stop("Invalid experiment")
    }
    exp[, grep(tag, exp$type, invert = TRUE, ignore.case = TRUE)]
}

#' @title Filter vector
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param vec Vector to be filtered
#' @param min Minimum value
#' @param max Maximum value
#' @param include.na Should NA's be included?
#' @noRd
filterVec <- function(vec, min, max, include.na = FALSE) {
    which(as.vector(vec >= min & vec < max) | is.infinite(vec) * include.na |
        is.na(vec) * include.na | is.nan(vec) * include.na)
}

#' @title Filter the compounds based on batch-corrected RSDQC
#' @description This function can be used to only keep compounds with
#' a RSDQC threshold. This ensures that only reliably measured compounds
#' are kept.
#' @details placeholder
#' @returns SummarizedExperiment with removed compounds based on the RSDQC
#' threshold provided
#' @param exp SummarizedExperiment object
#' @param min Minimum RSDQC value, defaults to 0
#' @param max Maximum RSDQC value, defaults to 30
#' @param include.na Should NA's be included?
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
#' # Perform analysis
#' exp <- doAnalysis(exp)
#'
#' # Keep compounds with a maximum RSDQC of 30
#' filterRSDQC(exp, min = 0, max = 30)
#'
#' # Keep compounds between RSDQC values 5 and 15
#' filterRSDQC(exp, min = 5, max = 15)
filterRSDQC <- function(exp, min = 0, max = 30, include.na = FALSE) {
    rsd <- rowData(exp)$rsdqcCorrected
    exp[filterVec(rsd, min, max, include.na), ]
}

#' @title Filter the Background Signal
#' @description This function can be used to only keep compounds with a
#' low background signal.
#' @details placeholder
#' @returns SummarizedExperiment with removed compounds based on the background
#' signal values and threshold provided
#' @param exp SummarizedExperiment object
#' @param min Minimum background signal value, defaults to 0
#' @param max Maximum background signal value, defaults to 0.4
#' @param type Aliquot type to use for background signal calculation. Defaults
#' to "BLANK"
#' @param include.na Should NA's be included?
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
#' # Filter by background signal
#' filterBackground(exp, min = 0, max = 0.4)
filterBackground <- function(exp, min = 0, max = 0.4, type = "BLANK",
                             include.na = FALSE) {
    be <- rowData(exp)$backgroundSignal
    be[is.nan(be)] <- NA
    exp[filterVec(be, min, max, include.na), ]
}

#' @title Filter aliquot QC outliers from the dataset
#' @description This function identifies and removes outlier QC aliquots from
#' the Experiment. Outliers are found by applying the RosnerTest on the median
#' A/IS ratio of QCs.
#' @description placeholder
#' @returns SummarizedExperiment with removed QC aliquots based their median
#' values
#' @param exp SummarizedExperiment object
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
#' # Remove QC outliers
#' filterOutliers(exp)
filterOutliers <- function(exp) {
    exp <- identifyOutliers(exp)

    return(exp[, exp$use])
}
