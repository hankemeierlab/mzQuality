#' @title Construct a Common Dataframe from Vendor File(s)
#' @description This function processes vendor files, which may come in
#'   different export formats, and constructs a common, human-readable
#'   table format. It ensures that all column names are standardized,
#'   replaces any value equal to 'QC' with 'SQC', and adds a batch number
#'   if none is present. Regular expressions are used to parse aliquot names.
#' @details The function handles vendor-specific formats by applying
#'   appropriate parsing routines. For Sciex files, it uses a specific
#'   routine to process the data. Missing values are replaced with `NA`,
#'   and column names are converted to lowercase for consistency.
#' @returns A data.frame in a generic, common long format with standardized
#'   column names and additional metadata such as batch numbers.
#' @param files A character vector of paths to batch files or a combined file.
#' @param vendor An optional character string specifying the vendor used.
#' @param regex A regular expression for parsing aliquot names.
#' @importFrom utils read.delim
#' @export
#' @examples
#' # Read example dataset
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
readData <- function(files, vendor = NA, regex = NULL) {
    mandatoryColumns <- c(
        "Acquisition Date & Time", "Acquisition.Date...Time",
        "Component Name", "Component.Name",
        "Retention Time", "Retention.Time",
        "Area"
    )

    df <- .setBatches(lapply(files, function(file) {
        df <- read.delim(file, sep = "\t")

        if (sum(colnames(df) %in% mandatoryColumns) > 1) {
            # Sciex routine
            df <- .processSciex(df, regex = regex)
        }

        df$datetime <- as.character(df$datetime)

        df[df == "N/A"] <- NA
        colnames(df) <- tolower(colnames(df))
        df$type <- toupper(df$type)
        df
    }))

    df
}

#' @title Get columns that contain metadata in a long dataframe format
#' @description Identifies columns in a dataframe that contain metadata based
#'     on unique combinations of values in the specified index column.
#' @details This function checks each column in the dataframe (excluding the
#'     index column) to determine if it contains metadata. A column is
#'     considered metadata if it has unique values for each unique value in
#'     the index column.
#' @param df A dataframe containing only columns with more than one unique
#'     value.
#' @param index A string specifying the column to use as the index.
#' @return A character vector of column names that contain metadata.
#' @importFrom dplyr distinct
#' @noRd
.getMetaData <- function(df, index) {
    cols <- c(index)
    for (col in colnames(df)) {
        if (col != index) {
            uniq_comb <- dplyr::distinct(df[, c(index, col)])
            if (!any(duplicated(uniq_comb[, index]))) {
                cols <- c(cols, col)
            }
        }
    }
    return(cols)
}

#' @title Set the metadata of compound or aliquot data
#' @description Creates or updates a metadata DataFrame for compounds or
#'     aliquots, ensuring proper row ordering and handling of missing rows.
#' @details This function ensures that the metadata DataFrame is properly
#'     formatted and ordered based on the provided rows. If the metadata
#'     DataFrame is empty, it initializes it with the unique rows. Otherwise,
#'     it updates the rownames and removes the specified index column.
#' @param rowData A DataFrame containing the metadata.
#' @param rowIndex A string specifying the column to use as the row index.
#' @param rows A character vector of row names to include in the metadata.
#' @return A DataFrame containing the updated metadata, ordered by rownames.
#' @importFrom dplyr distinct
#' @importFrom S4Vectors DataFrame
#' @noRd
.setMetaDataFrame <- function(rowData, rowIndex, rows) {
    if (nrow(rowData) == 0) {
        rowData <- DataFrame(row.names = unique(rows))
    } else {
        rownames(rowData) <- rows
        rowData <- rowData[-which(colnames(rowData) == rowIndex)]
    }
    return(rowData)
}

#' @title Build metadata DataFrame using an index column
#' @description Constructs a metadata DataFrame from a long-format dataframe
#'     using the specified index column.
#' @details This function extracts metadata columns from the input dataframe
#'     using `.getMetaData()`, then formats and orders the metadata using
#'     `.setMetaDataFrame()`.
#' @param df A dataframe in long format containing the data.
#' @param index A string specifying the column to use as the index.
#' @return A DataFrame containing the metadata, ordered by the unique values
#'     in the index column.
#' @importFrom dplyr distinct
#' @noRd
.buildMetadataWithIndex <- function(df, index) {
    rows <- unique(df[, index])
    rowData <- dplyr::distinct(df[.getMetaData(df, index)])
    rowData <- .setMetaDataFrame(rowData, index, rows)
    return(rowData)
}

#' @title Build assay list from a long-format dataframe
#' @description Constructs a list of assay matrices from a long-format
#'     dataframe, ensuring proper handling of missing combinations.
#' @details This function creates assay matrices for each assay column in the
#'     dataframe. It ensures that all combinations of rows and columns are
#'     present, filling in missing values if necessary. The resulting assay
#'     matrices are named after the assay columns.
#' @param df A dataframe in long format containing the data.
#' @param rowData A DataFrame containing row metadata.
#' @param colData A DataFrame containing column metadata.
#' @param rowIndex A string specifying the column to use as the row index.
#' @param colIndex A string specifying the column to use as the column index.
#' @return A named list of assay matrices, with rows corresponding to
#'     `rowData` and columns corresponding to `colData`.
#' @noRd
.buildAssayList <- function(df, rowData, colData, rowIndex, colIndex) {
    hasMissing <- nrow(df) != nrow(colData) * nrow(rowData)
    reserved <- c(colnames(rowData), colnames(colData), rowIndex, colIndex)
    assayNames <- setdiff(colnames(df), reserved)

    if (hasMissing) {
        df <- .fillMissingCombinations(
            df, rowIndex, colIndex, rowData, colData, assayNames
        )
    }

    assays <- lapply(assayNames, function(assayName) {
        m <- NA
        values <- df[, assayName, drop = TRUE]
        if (all(is.numeric(values))) {
            m <- matrix(
                data = as.double(values),
                ncol = nrow(colData),
                nrow = nrow(rowData),
                dimnames = list(rownames(rowData), rownames(colData)),
                byrow = TRUE
            )
        }

        return(m)
    })
    names(assays) <- assayNames
    assays <- assays[vapply(assays, length, numeric(1)) > 1]
    return(assays)
}


#' @title Create a SummarizedExperiment Object from a Long Dataframe
#' @description mzQuality works with a SummarizedExperiment object to perform
#'   its analysis. This function aids in that process by constructing one for
#'   you. It uses cardinality to identify in which slot a column should be
#'   stored. After creation, it will validate if the object is suitable for
#'   mzQuality to use.
#' @param df A dataframe in a long format. Should contain at least a column for
#'   compounds, aliquots, area, type, and datetime. See the vignette for an
#'   example of a valid dataframe.
#' @param rowIndex The column name in `df` that should be used for compounds.
#'   Defaults to 'compound'.
#' @param colIndex The column name in `df` that should be used for aliquots.
#'   Defaults to 'aliquot'.
#' @param primaryAssay The column name in `df` that represents the area or
#'   intensities of measured compounds. Defaults to 'area'.
#' @param secondaryAssay The column name in `df` that represents the area or
#'   intensities of measured internal standards. If no internal standards are
#'   present, this parameter can be omitted, and the A/IS ratio will
#'   be set to 1.
#' @param typeColumn The column name in `df` that specifies the type of each
#'   sample (e.g., 'SAMPLE', 'QC'). Defaults to 'type'.
#' @param qc The name of the QC sample used for batch correction found in the
#'   'type' column of `df`. Defaults to 'SQC'.
#' @param secondaryIndex The column name in `df` that should be used for
#'   internal standard compounds. Defaults to 'compound_is'.
#' @importFrom dplyr distinct
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames<-
#'   assays<- assay
#' @importFrom lubridate now
#' @returns A SummarizedExperiment object with dimensions compounds x aliquots.
#' @export
#' @examples
#' # Read example dataset
#' data <- readData(system.file(package = "mzQuality", "example.tsv"))
#'
#' # Construct experiment
#' exp <- buildExperiment(
#'     df = data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area_is"
#' )
#'
#' # Construct experiment without concentrations
#' data <- data[-which(colnames(data) == "concentration")]
#' exp <- buildExperiment(
#'     df = data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area_is"
#' )
#' isValidExperiment(exp)
#'
#' # Construct experiment without internal standards
#' data <- data[-grep("_is", colnames(data))]
#' exp <- buildExperiment(
#'     df = data,
#'     rowIndex = "compound",
#'     colIndex = "aliquot",
#'     primaryAssay = "area",
#'     secondaryAssay = "area"
#' )
buildExperiment <- function(
        df, rowIndex = "compound", colIndex = "aliquot",
        primaryAssay = "area", secondaryAssay = "area_is",
        typeColumn = "type", qc = "SQC", secondaryIndex = "compound_is"
) {
    df <- as.data.frame(df)
    colnames(df) <- tolower(colnames(df))

    if (!"type" %in% colnames(df)) df$type <- "SAMPLE"
    stopifnot(isValidDataframe(df))

    df$type <- toupper(df$type)

    if (!qc %in% df$type) {
        message(sprintf("QC '%s' not found, default to type 'SAMPLE'", qc))
        qc <- "SAMPLE"
    }

    if (!secondaryAssay %in% colnames(df)) secondaryAssay <- primaryAssay

    allNa <- colSums(is.na(df)) == nrow(df)
    allSame <- vapply(df, function(col) all(col == col[1]), logical(1))
    keep <- !allNa & (!allSame | is.na(allSame))

    df <- df[, keep, drop = FALSE]

    idx <- colnames(df) != typeColumn
    rowData <- .buildMetadataWithIndex(df[, idx], rowIndex)
    colData <- .buildMetadataWithIndex(df, colIndex)
    assays <- .buildAssayList(df, rowData, colData, rowIndex, colIndex)

    hasIS <- secondaryAssay %in% names(assays) & primaryAssay != secondaryAssay
    meta <- list(
        QC = qc,
        primary = primaryAssay,
        secondary = secondaryAssay,
        hasIS = hasIS,
        Date = Sys.time()
    )

    exp <- SummarizedExperiment(
        assays = assays, rowData = rowData, colData = colData,
        metadata = meta
    ) %>%
        .addConcentrationMetadata() %>%
        .addInitialAnalysis()

    return(exp)
}

#' @title Set the Colors of the Samples
#' @description Assigns unique colors to sample types using the viridis
#'   color palette. Ensures that each sample type is assigned a distinct
#'   color for visualization purposes.
#' @details This function uses the viridis color palette to generate a
#'   range of colors for the provided sample types. It ensures that the
#'   colors are consistent and can be overwritten if additional arguments
#'   are provided.
#' @returns A named character vector of colors corresponding to the sample
#'   types.
#' @param types A character vector with the types of the samples.
#' @param ... Additional arguments to be passed to the function.
#' @importFrom viridis viridis
#' @importFrom stats setNames
#' @noRd
.setSampleColors <- function(types, ...) {
    types <- unique(types)
    colorVec <- viridis(length(types), begin = 0.2, end = 0.8, option = "H")

    colorVec <- setNames(colorVec, types)
    colorVec <- colorVec[setdiff(types, names(.sampleColors()))]

    colorVec <- .sampleColors(colorVec, ..., overwrite = length(list(...)) > 0)
    return(colorVec)
}


#' @title Finalize a SummarizedExperiment Object
#' @description Finalizes the creation of a SummarizedExperiment object by
#'   adding metadata, ordering samples, and performing additional
#'   calculations.
#' @details This function orders the samples by injection time, assigns
#'   colors to sample types, and calculates additional metadata such as
#'   batch and order. It also performs various calculations, including
#'   identifying outliers and adding linear ranges.
#' @returns A finalized SummarizedExperiment object with additional metadata
#'   and calculated values.
#' @param exp A SummarizedExperiment object.
#' @importFrom viridis viridis
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom lubridate format_ISO8601
#' @importFrom dplyr group_by n pull arrange mutate
#' @noRd
.addInitialAnalysis <- function(exp) {
    exp <- exp[, order(exp$injection_time, method = "radix")]
    exp$datetime <- lubridate::as_datetime(exp$injection_time)
    exp$color <- .setSampleColors(exp$type)[exp$type]

    batches <- rep(1, ncol(exp))
    if ("batch" %in% colnames(colData(exp))) {
        batches <- as.integer(as.factor(exp$batch))
    }

    batches <- formatC(batches, max(nchar(batches)), format = "d", flag = "0")
    exp$batch <- paste0("Batch", batches)

    if (!"order" %in% colnames(colData(exp))) {
        order <- colData(exp) %>%
            as.data.frame() %>%
            arrange(.data$injection_time) %>%
            group_by(.data$batch) %>%
            mutate(order = seq_len(n())) %>%
            pull(order)

        colData(exp)$order <- order
    }

    exp <- calculateRatio(exp) %>%
        identifyOutliers() %>%
        identifyMisInjections() %>%
        addLinearRange(calType = "CAL", saveAssay = "CalRange") %>%
        addLinearRange(calType = "ACAL", saveAssay = "ACALRange")

    return(exp)
}

#' @title Fill Missing Aliquot-Compound Combinations
#' @description Identifies and fills in missing combinations of aliquots
#'   and compounds in a long-format dataframe.
#' @details This function ensures that all possible combinations of aliquots
#'   and compounds are present in the dataframe. Missing combinations are
#'   added with NA values for the assay columns.
#' @returns A dataframe with all aliquot-compound combinations, including
#'   the missing ones filled with NA values.
#' @param df A dataframe in long format.
#' @param rowIndex The column name in `df` to be used for compounds.
#' @param colIndex The column name in `df` to be used for aliquots.
#' @param rowData A dataframe containing row data (e.g., compounds).
#' @param colData A dataframe containing column data (e.g., aliquots).
#' @param assayNames A character vector of column names in `df` representing
#'   the area or intensities of measured compounds.
#' @importFrom dplyr bind_rows
.fillMissingCombinations <- function(df, rowIndex, colIndex, rowData, colData,
                                    assayNames) {

    # Identify compounds and aliquots with missing combinations
    incompleteCompounds <- names(which(table(df[[rowIndex]]) < nrow(colData)))
    incompleteAliquots <- names(which(table(df[[colIndex]]) < nrow(rowData)))

    # All possible missing pairs to add
    grid <- expand.grid(
        aliquot = incompleteAliquots,
        compound = incompleteCompounds
    )

    pairs <- paste0(grid$aliquot, grid$compound)

    # Existing compound-aliquot pairs
    existingPairs <- paste0(df[[colIndex]], df[[rowIndex]])

    # Keep only the missing ones
    grid <- grid[!pairs %in% existingPairs, ]
    colnames(grid) <- c(colIndex, rowIndex)

    # Fill assay columns with NA
    for (assay in assayNames) {
        grid[[assay]] <- NA
    }

    # Combine original and missing rows
    df <- bind_rows(df[, c(colIndex, rowIndex, assayNames)], grid)

    # Sort for consistency
    df <- df[order(df[[colIndex]], df[[rowIndex]], method = "radix"), ]

    return(df)
}

#' @title Convert a SummarizedExperiment to a data frame in the long format
#' @description A SummarizedExperiment is great for analyzing data, but not
#' useful for storage in human-readable format. To facilitate in this need,
#' this function can be used to convert a SummarizedExperiment to a dataframe
#' in the long format. Note that metadata won't be exported.
#' @param exp a SummarizedExperiment object
#' @param rowIndex Columnname in the resulting dataframe for compounds,
#' defaults to 'Compound'.
#' @param colIndex Columnname in the resulting dataframe for aliquots,
#' defaults to 'Aliquot'.
#' @importFrom SummarizedExperiment assayNames assay colData rowData
#' @importFrom utils stack
#' @returns a data.frame in a generic, human-readable format
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
#' # Convert to a table format
#' expToCombined(exp, rowIndex = "Compound", colIndex = "Aliquot")
expToCombined <- function(exp,  rowIndex = "compound", colIndex = "aliquot") {
        # transform to a long format, only maintaining the dim-names

        rData <- rowData(exp)
        cData <- colData(exp)

        x <- expand.grid(list(rownames(rData), rownames(cData)))
        colnames(x) <- c(rowIndex, colIndex)

        # Get long-format dataframe of all assays in the experiment
        assays <- do.call(cbind, lapply(assayNames(exp), function(assay){
            # Per assay, retrieve it as a matrix
            mat <- assay(exp, assay)
            return(stack(as.data.frame(mat))[, -2, drop = FALSE])
        }))

        # Set the column names of the assays to their assay names
        colnames(assays) <- assayNames(exp)

        # Set a vector of all columns names to be used after binding
        cols <- c(
            colnames(x),
            colnames(assays),
            colnames(cData),
            colnames(rData)
        )

        df <- as.data.frame(cbind(
            x,
            assays,
            cData[x[, colIndex], ],
            rData[x[, rowIndex], ]
        ))

        # Order the dataframe based on the Compounds and aliquots
        df <- df[order(df[, rowIndex], df[, colIndex]), , drop = FALSE]

        # Assign new rownames
        rownames(df) <- seq_len(nrow(df))

        # Assign the correct column names
        colnames(df) <- cols

        # Return the data.fame
        cols <- !duplicated(colnames(df)) & !is.na(colnames(df))
        return(df[, cols, drop = FALSE])
}

#' @title Add Known Concentrations of Calibration Lines
#' @description Adds known concentrations of calibration lines (CAL) to a
#'     SummarizedExperiment object. This function creates a new assay called
#'     "Concentration" in the experiment, which is filled with `NA` values
#'     and the supplied concentrations.
#' @details The concentrations should be provided as a data frame where rows
#'     represent the names of calibration lines and columns represent
#'     compound names. Optionally, the `filterComps` argument can be set to
#'     TRUE to retain only compounds in the experiment that have a known
#'     concentration in the CAL lines.
#' @param exp A SummarizedExperiment object to which the concentration data
#'     will be added.
#' @param df A data frame containing the concentrations of calibration lines.
#'     Rows should correspond to calibration line names, and columns should
#'     correspond to compound names.
#' @returns A SummarizedExperiment object with an added assay called
#'     "Concentration". The assay is filled with the values provided in `df`,
#'     and other values are set to `NA`.
#' @export
#' @examples
#' # Read example dataset
#' exp <- readRDS(system.file(package = "mzQuality", "data.RDS"))
#' # Add concentrations
#' file <- system.file(package = "mzQuality", "concentrations.txt")
#' concentrations <- read.delim(file)
#'
#' addConcentrations(exp, df)
addConcentrations <- function(exp, df) {
    comps <- intersect(rownames(exp), rownames(df))

    find <- gregexec("([A-Za-z]+)([0-9]+)", colnames(df))
    matches <- do.call(rbind, lapply(regmatches(colnames(df), find), t))
    idx1 <- exp$calno %in% as.integer(matches[, 3])
    idx2 <- exp$type %in% matches[, 2]
    aliqs <- colnames(exp)[idx1 & idx2]

    if (length(comps) == 0 || length(aliqs) == 0) return(exp)

    m <- assay(exp)
    m[TRUE] <- NA

    # Subset concentration when not all expected callines are present
    cols <- aliqs[aliqs %in% colnames(exp[, exp$batch == unique(exp$batch)[1]])]
    df <- df[, exp[, cols]$calno + 1]

    res <- do.call(cbind, lapply(unique(exp$batch), function(batch) {
        cols <- intersect(aliqs, colnames(exp[, exp$batch == batch]))
        colnames(df) <- cols
        df
    }))

    m[comps, aliqs] <- as.matrix(res[comps, aliqs])
    columns <- exp$type == metadata(exp)$concentration & exp$calno == 0

    if (sum(columns) > 0) {
        concs <- m[, columns]
        concs[is.na(concs)] <- 0
        m[, columns] <- concs
    }

    assay(exp, "concentration") <- m
    exp <- .addConcentrationMetadata(exp)
    return(exp)
}

#' @title Add metadata when concentrations are supplied
#' @description Adds metadata to a SummarizedExperiment object when
#' concentrations are supplied. This function checks if the concentrations are
#' present in the assay names and adds metadata accordingly.
#' @details The function checks if the "concentration" assay is present in the
#' experiment. If it is, it retrieves the concentration type and updates the
#' metadata and rowData of the experiment. It also sets flags for known and
#' calculated concentrations.
#' @param exp A SummarizedExperiment object to which the concentration
#' metadata will be added.
#' @returns A SummarizedExperiment object with updated metadata and rowData.
#' @importFrom SummarizedExperiment assayNames assay rowData rowData<-
#' @noRd
.addConcentrationMetadata <- function(exp) {
    if (!"concentration" %in% assayNames(exp)) {
        return(exp)
    }

    m <- assay(exp, "concentration")
    concType <- unique(exp$type[colSums(is.na(m)) != nrow(m)])
    if (length(concType) > 1) {
        message("> 1 calibration line is not supported, selecting first type")
        concType <- concType[1]
    }

    metadata(exp)$concentration <- concType
    rowData(exp)$hasKnownConcentrations <- FALSE
    hasConcentrations <- rowSums(m, na.rm = TRUE) > 0
    rowData(exp[hasConcentrations, ])$hasKnownConcentrations <- TRUE

    return(exp)
}

#' @title Calculate the A/IS ratio with minimal imputation
#' @param exp SummarizedExperiment object
#' @param assay1 Name of the primary assay, usually "Area" for compounds
#' @param assay2 Name of the secondary assay, usually "Area_is" for internal
#' standards
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay assay<-
#' @returns SummarizedExperiment with an assay slot called "ratio" containing
#' the ratio between two assays.
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
#' # Calculate the ratio
#' exp <- calculateRatio(exp)
#'
#' # Show results
#' assay(exp, "ratio")
calculateRatio <- function(
        exp, assay1 = metadata(exp)$primary,
        assay2 = metadata(exp)$secondary
) {
    if (!metadata(exp)$hasIS) {
        df <- assay(exp, assay1)
    } else {
        df <- assay(exp, assay1) / assay(exp, assay2)
    }

    df[df <= 0] <- NA
    assay(exp, "ratio") <- df

    return(exp)
}
