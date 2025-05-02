#' @title Validate Measurements
#' @description Validates if a dataframe is suitable for use with the
#'   `mzQuality` package.
#' @details This function checks if the dataframe contains the required
#'   columns, has no duplicate aliquot-batch-compound combinations, and
#'   includes at least one "QC" type.
#' @returns A boolean value indicating whether the dataframe is valid.
#' @param dataframe A data.frame containing measurements.
#' @importFrom dplyr %>%
#' @export
#' @examples
#' # Read example dataset
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
#'
#' # Validate the dataframe
#' isValid <- isValidDataframe(data)
isValidDataframe <- function(dataframe = NULL) {
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

    x <- data.frame(
        c = make.names(dataframe$compound),
        a = make.names(dataframe$aliquot),
        b = make.names(dataframe$batch)
    )

    if (any(duplicated(x))) {
        x <- x[which(duplicated(x)), ]
        message("Duplicated aliquot-batch-compound combination!")
        return(FALSE)
    }

    return(TRUE)
}

#' @title Validate Experiment Object
#' @description Validates if a SummarizedExperiment object is correctly
#'   formatted and suitable for use with the `mzQuality` package.
#' @details This function checks that the object is a valid
#'   SummarizedExperiment, contains the required metadata columns in
#'   `colData`, and includes the "ratio" assay.
#' @returns A boolean value indicating whether the SummarizedExperiment
#'   object is valid.
#' @param exp A SummarizedExperiment object.
#' @importFrom methods is
#' @export
#' @examples
#' # Read example dataset
#' exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
#'
#' # Validate the experiment
#' isValid <- isValidExperiment(exp)
isValidExperiment <- function(exp) {
    if (!is(exp, "SummarizedExperiment")) return(FALSE)

    valid <- all(
        nrow(exp) > 0,
        ncol(exp) > 0,
        all(c("type", "datetime", "batch") %in% colnames(colData(exp))),
        "ratio" %in% assayNames(exp)
    )

    return(valid)
}

#' @title Convert an experiment to one capable with mzQuality
#' @description This function can be used to convert an existing
#' SummarizedExperiment to one capable to use with mzQuality. This function
#' will change some colData and rowData columns, but will leave assays
#' untouched.
#' @details This function will convert the column names of the colData and
#' rowData to the ones used in mzQuality. The function will also add
#' metadata to the SummarizedExperiment object, which is used to store
#' information about the experiment. The metadata will be added to the
#' SummarizedExperiment object if it is not already present. The function
#' will also check if the primaryAssay and secondaryAssay are the same.
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
#' @importFrom SummarizedExperiment colData<-
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
#' convertExperiment(exp, primaryAssay = "Area")
convertExperiment <- function(
    exp, primaryAssay = "area", secondaryAssay = primaryAssay,
    istd = "compound_is", type = "type", datetime = "datetime",
    batch = "batch", qcType = "SQC"
) {

    cols <- colnames(rowData(exp))

    if (primaryAssay != secondaryAssay) {
        cols[which(cols) == istd] <- "compound_is"
        colnames(rowData(exp)) <- cols
    }

    cols <- colnames(colData(exp))
    cols[which(cols == type)] <- "type"
    cols[which(cols == datetime)] <- "datetime"
    cols[which(cols == batch)] <- "batch"

    x <- colData(exp)
    colnames(x) <- cols
    colData(exp) <- x

    mandatory <- list(
        QC = qcType, primary = primaryAssay,
        secondary = secondaryAssay, hasIS = primaryAssay != secondaryAssay,
        Date = Sys.time()
    )

    toAdd <- setdiff(names(mandatory), names(metadata(exp)))
    if (length(toAdd) > 0) {
        metadata(exp)[toAdd] <- mandatory[toAdd]
    }

    .addInitialAnalysis(exp)
}
