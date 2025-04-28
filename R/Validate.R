#' @title Validate measurements
#' @description placeholder
#' @details placeholder
#' @returns Boolean value if the dataframe is valid to use with mzQuality
#' @param dataframe data.frame with measurements.
#' @importFrom dplyr %>%
#' @export
#' @examples
#' # Read example dataset
#' data <- read.delim(system.file(package = "mzQuality", "dataset.txt"))
#'
#' # Validate if dataframe is corrct
#' validateDataframe(data)
validateDataframe <- function(dataframe = NULL) {
    columns <- c("aliquot", "compound", "type", "batch", "datetime")
    check <- all(
        !is.null(dataframe),
        nrow(dataframe) > 0,
        ncol(dataframe) > 0,
        columns %in% colnames(dataframe)
    )

    if (!check) {
        return(FALSE)
    }

    uts <- unique(dataframe$type)
    if (sum(grepl("QC", uts)) == 0) {
        return(FALSE)
    }

    x <- data.frame(
        c = make.names(dataframe$compound),
        a = make.names(dataframe$aliquot),
        b = make.names(dataframe$batch)
    )

    if (any(duplicated(x))) {
        x <- x[which(duplicated(x)), ]
        stop("Duplicated aliquot-batch-compound!")
    }

    return(TRUE)
}

#' @title Validate experiment object
#' @description placeholder
#' @details placeholder
#' @returns Boolean value indicating if the SummarizedExperiment is valid
#' to use with mzQuality
#' @param exp SummarizedExperiment object
#' @export
#' @importFrom methods is
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
#' # Validate if the experiment is formatted correctly
#' validateExperiment(exp)
validateExperiment <- function(exp) {
    all(
        methods::is(exp, "SummarizedExperiment"),
        nrow(exp) > 0,
        ncol(exp) > 0,
        all(c("type", "datetime", "batch") %in% colnames(colData(exp))),
        "ratio" %in% assayNames(exp)
    )
}

#' @title Validate calibration dataframe
#' @description TO ADD
#' @param exp SummarizedExperiment object
#' @param df Dataframe with concentration measured in calibration lines
#' @returns Boolean value indicating if the calibration dataframe is suitable
#' for use with mzQuality
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
#' # Add calibration file + checks
validateCalibration <- function(exp, df) {
    # Validate that it is correct
    find <- gregexec("([a-zA-Z]+)([0-9]+)", colnames(df))
    matches <- regmatches(colnames(df), find)
    matches <- do.call(rbind, lapply(matches, t))
    idx1 <- exp$calno %in% as.integer(matches[, 3])
    idx2 <- exp$type %in% matches[, 2]
    aliqs <- colnames(exp)[idx1 & idx2]
    all(
        any(aliqs %in% colnames(exp)),
        any(rownames(df) %in% rownames(exp))
    )
}

#' @title Convert an experiment to one capable with mzQuality
#' @description This function can be used to convert an existing
#' SummarizedExperiment to one capable to use with mzQuality. This function
#' will change some colData and rowData columns, but will leave assays
#' untouched.
#' @details
#' @returns SummarizedExperiment that is suitable for use with mzQuality
#' @param exp A SummarizedExperiment object
#' @param primaryAssay Assay that should be used for compound Areas / Intensity
#' @param secondaryAssay Assay that should be used for internal standard
#' Areas / Intensity. Defaults to equal to primaryAssay, meaning no internal
#' standards are present.
#' @param istd Column in the rowData that represents the internal standard used
#' for the compound. Defaults to 'compound_is'
#' @param type Column in the colData that represents the sample type.
#' Defaults to 'Type'
#' @param datetime Column in the colData that represents the datetime when the
#' sample was measured. Defaults to 'Datetime'
#' @param batch Column in the colData that represents the batch name / number.
#' Defaults to 'Batch'
#' @param qcType Name of the QC sample type used. Defaults to 'SQC'
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
#' convertExperiment(exp, primaryAssay = "Area")
convertExperiment <- function(
    exp, primaryAssay = "area", secondaryAssay = primaryAssay,
    istd = "compound_is", type = "type", datetime = "datetime",
    batch = "batch", qcType = "SQC") {
    if (requireNamespace("methods", quietly = TRUE)) {
        if (!methods::is(exp, "SummarizedExperiment")) {
            logwarn("Object is not of class 'SummarizedExperiment'")
            return(NULL)
        }
    }


    assays <- assayNames(exp)
    if (!primaryAssay %in% assays) {
        stop(sprintf("Cannot find assay '%s' in assayNames", primaryAssay))
    }
    if (!secondaryAssay %in% assays) {
        stop(sprintf("Cannot find assay '%s' in assayNames", secondaryAssay))
    }

    cols <- colnames(rowData(exp))
    if (!istd %in% cols) {
        stop(sprintf("Cannot find column '%s' in rowData", istd))
    }
    if (primaryAssay != secondaryAssay) {
        cols[which(cols) == istd] <- "compound_is"
        colnames(rowData(exp)) <- cols
    }

    cols <- colnames(colData(exp))

    if (!type %in% cols) {
        stop(sprintf("Cannot find column '%s' in colData", type))
    }
    if (!datetime %in% cols) {
        stop(sprintf("Cannot find column '%s' in colData", datetime))
    }

    cols[which(cols == type)] <- "type"
    cols[which(cols == datetime)] <- "datetime"
    cols[which(cols == batch)] <- "batch"

    if (!qcType %in% exp$type) {
        stop(sprintf(
            "Cannot find type '%s' in the column '%s'",
            qcType, type
        ))
    }

    x <- colData(exp)
    colnames(x) <- cols
    colData(exp) <- x

    metadata(exp) <- c(metadata(exp), list(
        QC = qcType, primary = primaryAssay,
        secondary = secondaryAssay, hasIS = primaryAssay != secondaryAssay,
        Date = lubridate::now()
    ))
    metadata(exp) <- metadata(exp)[duplicated(names(metadata(exp)))]
    finishExperiment(exp)
}
