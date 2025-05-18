#' @title Filter compounds that are equal to given internal standards
#' @description In a scenario where internal standards were measured as
#' compounds, this function aids in removing these from the
#' SummarizedExperiment before doing analysis. It filter by both the `tag`
#' argument and by checking if the A/IS ratio equals 1.
#' @details This function is useful when internal standards are
#' measured as compounds and you want to remove them from the analysis.
#' It filters the compounds based on the provided tag and checks if the A/IS
#' ratio equals 1. If the `hasIS` metadata is set to TRUE, it will also
#' remove compounds that have the same A/IS ratio as the internal standard.
#' @returns SummarizedExperiment with removed compounds based on the filter
#' provided
#' @param exp A SummarizedExperiment object
#' @param tag Tag to search for in compound names that indicates an internal
#' standard. Defaults to 'ISTD'
#' Regular expressions are supported.
#' @export
#' @examples
#' # Read example dataset
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
#'
#' # Construct experiment
#' exp <- buildExperiment(
#'     data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area_is"
#' )
#'
#' # Filter internal standards from compounds
#' filterISTD(exp, tag = "ISTD")
#'
#' # Filter all standards from compounds
#' filterISTD(exp, tag = "STD")
filterISTD <- function(exp, tag = "STD") {
    stopifnot(isValidExperiment(exp))

    exp <- exp[grep(tag, rownames(exp), invert = TRUE, ignore.case = TRUE), ]
    if (metadata(exp)$hasIS) {
        a <- assay(exp, metadata(exp)$primary)
        b <- assay(exp, metadata(exp)$secondary)
        numSame <- rowSums(a / b == 1 | is.na(a / b), na.rm = TRUE)
        exp <- exp[numSame != ncol(exp), ]
    }
    exp
}

#' @title Filter SST Aliquots from the SummarizedExperiment
#' @description System Sustainability Tests (SST) samples are often executed by
#' the technician for control during measurements. However, they are often
#' unwanted in analysis of measurements. This function helps to remove these
#' SSTs from the SummarizedExperiment object
#' @details This function is useful when you want to remove SST aliquots from
#' the analysis. It filters the aliquots based on the provided tag.
#' The tag is used to search for SST aliquots in the `type` column of the
#' SummarizedExperiment object. The function uses regular expressions to
#' match the tag, so you can use any pattern that matches the SST aliquots.
#' @returns SummarizedExperiment with removed aliquots based on the filter
#' provided
#' @param exp A SummarizedExperiment object
#' @param tag Tag to search for in aliquot types that indicate a SST.
#' Defaults to 'SST'
#' Regular expressions are supported.
#' @export
#' @examples
#' # Read example dataset
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
#'
#' # Construct experiment
#' exp <- buildExperiment(
#'     data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area_is"
#' )
#'
#' # Filter SST aliquots
#' filterSST(exp, tag = "SST")
filterSST <- function(exp, tag = "SST") {
    stopifnot(isValidExperiment(exp))
    exp[, grep(tag, exp$type, invert = TRUE, ignore.case = TRUE)]
}

#' @title Filter a Vector
#' @description Filters a vector based on minimum and maximum values, with an
#'   option to include NA values.
#' @details This function identifies indices of elements in the vector that
#'   fall within the specified range or meet the NA inclusion criteria.
#' @returns A vector of indices that meet the filtering criteria.
#' @param vec A vector to be filtered.
#' @param min The minimum value for filtering.
#' @param max The maximum value for filtering.
#' @param include.na Logical. Should NA values be included?
#' @noRd
.filterVec <- function(vec, min, max, include.na = FALSE) {
    which(as.vector(vec >= min & vec < max) | is.infinite(vec) * include.na |
        is.na(vec) * include.na | is.nan(vec) * include.na)
}

#' @title Filter the compounds based on batch-corrected RSDQC
#' @description This function can be used to only keep compounds with
#' a RSDQC threshold. This ensures that only reliably measured compounds
#' are kept.
#' @details This function is useful when you want to filter compounds based on
#' the RSDQC values. It filters the compounds based on the provided
#' minimum and maximum RSDQC values. The function uses regular expressions
#' to match the tag, so you can use any pattern that matches the RSDQC
#' values.
#' @returns SummarizedExperiment with removed compounds based on the RSDQC
#' threshold provided
#' @param exp SummarizedExperiment object
#' @param min Minimum RSDQC value, defaults to 0
#' @param max Maximum RSDQC value, defaults to 30
#' @param include.na Should NA's be included?
#' @export
#' @examples
#' # Read example dataset
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
#'
#' # Construct experiment
#' exp <- buildExperiment(
#'     data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area_is"
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
    exp[.filterVec(rsd, min, max, include.na), ]
}

#' @title Filter the Background Signal
#' @description This function can be used to only keep compounds with a
#' low background signal.
#' @details This function is useful when you want to filter compounds based on
#' the background signal values. It filters the compounds based on the
#' provided minimum and maximum background signal values.
#' @returns SummarizedExperiment with removed compounds based on the background
#' signal values and threshold provided
#' @param exp SummarizedExperiment object
#' @param min Minimum background signal value, defaults to 0
#' @param max Maximum background signal value, defaults to 0.4
#' @param include.na Should NA's be included?
#' @export
#' @examples
#' # Read example dataset
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
#'
#' # Construct experiment
#' exp <- buildExperiment(
#'     data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area_is"
#' )
#'
#' # Filter by background signal
#' filterBackground(exp, min = 0, max = 0.4)
filterBackground <- function(exp, min = 0, max = 0.4, include.na = FALSE) {
    be <- rowData(exp)$backgroundSignal
    be[is.nan(be)] <- NA
    exp[.filterVec(be, min, max, include.na), ]
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
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
#'
#' # Construct experiment
#' exp <- buildExperiment(
#'     data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area_is"
#' )
#'
#' # Remove QC outliers
#' filterOutliers(exp)
filterOutliers <- function(exp) {
    exp <- identifyOutliers(exp)

    return(exp[, exp$use])
}
