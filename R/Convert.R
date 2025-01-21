#' @title Construct a common Dataframe from vendor file(s)
#' @description Vendor files come in different export formats. This function
#' constructs a common, human-readable table format. It ensures that all
#' columnnames are titled, replaces any value equal to 'QC' to 'SQC' and adds
#' a batch number if none is present. Regex are used to parse aliquot names.
#' @details placeholder
#' @returns placeholder
#' @param files Path to batch files or combined file
#' @param vendor optional Character, which vendor is used?
#' @param regex Regular Expression for parsing the aliquots
#' @returns data.frame in a generic, common long format
#' @importFrom arrow read_delim_arrow
#' @export
buildCombined <- function(files, vendor = NA, regex = NULL) {
    if (length(files) == 1) {
        df <- arrow::read_delim_arrow(files[1], delim = "\t")

        mandatoryColumns <- c("Acquisition Date & Time", "Acquisition.Date...Time",
                              "Component Name", "Component.Name",
                              "Retention Time", "Retention.Time",
                              "Area"
        )

        if (sum(colnames(df) %in% mandatoryColumns) > 1) {
            df <- sciex(df, regex = regex)
            df[df == "N/A"] <- NA
            df <- setBatches(list(df))
        }

        colnames(df) <- tolower(colnames(df))

        df$type <- toupper(df$type)

    } else {
        df <- setBatches(lapply(files, function(file) {
            combined <- arrow::read_delim_arrow(file, delim = "\t")

            if (startsWith(colnames(combined)[1], "Sample")) {
                combined <- sciex(combined, regex = regex)
            }
            combined[combined == "N/A"] <- NA

            colnames(combined) <- tolower(colnames(combined))
            combined$type <- toupper(combined$type)
            combined
        }))
    }

    df
}

#' @title Get columns that contain metadata in a long dataframe format
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param df Dataframe with only columns with > 1 unique values
#' @param index Which index should be used as the index column
#' @importFrom dplyr distinct
getMetaData <- function(df, index) {
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

#' @title Set the Metadata of compound or aliquot data
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param rowData
#' @param rowIndex
#' @param rows
#' @importFrom dplyr distinct
setMetaDataFrame <- function(rowData, rowIndex, rows) {
  if (nrow(rowData) == 0) {
    rowData <- DataFrame(row.names = unique(rows))
  } else {
    rownames(rowData) <- rows
    rowData <- rowData[-which(colnames(rowData) == rowIndex)]
  }
  #return(rowData)#[order(rownames(rowData), method = "radix"), , drop = FALSE])
  return(rowData[order(rownames(rowData), method = "radix"), , drop = FALSE])
}

#' @title Create a SummarizedExperiment object from a long dataframe
#' @description mzQuality works with a SummarizedExperiment object to perform
#' its analysis. This function aids in that process by constructing one for
#' you. It uses cardinality to identify in which slot a column should be
#' stored. After creation, it will validate if the object is suitable for
#' mzQuality to use.
#' @param df Dataframe in a long format. Should contain at least a column for
#' compounds, aliquots, area, type and datetime. See the vignette for an
#' example of a valid dataframe.
#' @param rowIndex Column name in `df` that should be used for compounds.
#' Defaults to 'Compound'
#' @param colIndex Column name in `df` that should be used for aliquots.
#' Defaults to 'Aliquot'
#' @param primaryAssay Column name in `df` that represents the
#' area / intensities of measured compounds. Defaults to 'Area'
#' @param secondaryAssay Column name in `df` that represents the
#' area / intensities of measured internal standards.
#' If no internal standards are present, this parameter can be omitted and the
#'  A/IS ratio will be set to 1.
#' @param qc Name of the QC sample used for batch correction found in the
#' 'Type' column of `df`. Defaults to 'SQC'
#' @importFrom stringr str_to_title
#' @importFrom dplyr distinct
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames<- assays<- assay
#' @importFrom lubridate now
#' @returns A SummarizedExperiment object with dimensions compounds x aliquots.
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
#' exp
buildExperiment <- function(df,
                            rowIndex = "compound",
                            colIndex = "aliquot",
                            primaryAssay = "area",
                            secondaryAssay = "area_is",
                            typeColumn = "type",
                            qc = "SQC",
                            secondaryIndex = "compound_is") {
  df <- as.data.frame(df)

  colnames(df) <- tolower(colnames(df))
  df$type <- toupper(df$type)

  if (!qc %in% df$type) {

    if (length(unique(df$type)) == 1) {
      qc <- df$type[1]
    } else {
      stop("Could not find the given QC type")
    }
  }

  if (!any(c(rowIndex, colIndex) %in% colnames(df))) {
      stop(glue::glue("Cannot find columns {rowIndex} and {colIndex} in `df`",
                      rowIndex = rowIndex, colIndex = colIndex)
      )
  }

  if (!"type" %in% colnames(df)){
      stop("Cannot find mandatory column type in `df`")
  }


  if (!secondaryAssay %in% colnames(df)) {
    secondaryAssay <- primaryAssay
  } else {
     # If only one internal standard is used
  }


  rows <- unique(df[, rowIndex])
  cols <- unique(df[, colIndex])

  reserved <- df[, typeColumn]

  columns <- sapply(colnames(df), function(x){all(df[,x] == df[1,x]) | all(is.na(df[,x]))})
  columns[is.na(columns)] <- FALSE
  df <- df[, !columns]


  df$type <- reserved
  rowData <- dplyr::distinct(df[getMetaData(df, rowIndex)])
  colData <- dplyr::distinct(df[getMetaData(df, colIndex)])




  rowData <- rowData[, which(!colnames(rowData) %in% typeColumn), drop = FALSE]


  assays <- colnames(df)[!colnames(df) %in% c(
    colnames(rowData),
    colnames(colData),
    rowIndex,
    colIndex
  )]

  rowData <- setMetaDataFrame(rowData, rowIndex, rows)
  colData <- setMetaDataFrame(colData, colIndex, cols)

  df <- df[order(df[, colIndex], df[, rowIndex], method = "radix"), ]
  hasIS <- secondaryAssay %in% colnames(df) & primaryAssay != secondaryAssay

  if (hasIS && !secondaryIndex %in% colnames(rowData)) {
      # error, multiple internal standards per compound
  }

  hasMissing <- nrow(df) != nrow(colData) * nrow(rowData)
  message("Missing combinations found, trying to fix with NAs..")


  if (hasMissing) {
    df <- fillMissingCombinations(
      df, rowIndex, colIndex,
      rowData, colData, assays
    )
  }

  assays <- lapply(setNames(assays, assays), function(assayName) {
    matrix(data = as.double(df[, assayName]),
           ncol = nrow(colData),
           nrow = nrow(rowData),
           dimnames = list(rownames(rowData), rownames(colData))
    )
  })

  exp <- SummarizedExperiment(
    assays = assays,
    rowData = rowData,
    colData = colData,
    metadata = list(
      QC = qc,
      primary = primaryAssay,
      secondary = secondaryAssay,
      hasIS = hasIS,
      Date = lubridate::now()
    ))

  goodAssays <- vapply(assayNames(exp), function(x){
    all(is.numeric(assay(exp, x)))
  }, logical(1))

  assays(exp) <- lapply(names(which(goodAssays)), function(x) assay(exp, x))
  assayNames(exp) <- names(which(goodAssays))
  rownames(exp) <- make.names(rownames(exp))
  dimnames(exp) <- list(rownames(rowData), rownames(colData))


  finishExperiment(exp)
}

setSampleColors <- function(types, ...) {
    types <- unique(types)
    colorVec <- viridis::viridis(length(types), begin = 0.2, end = 0.8, option = "H")

    colorVec <- setNames(colorVec, types)
    colorVec <- colorVec[setdiff(types, names(sampleColors()))]

    colorVec <- sampleColors(colorVec, ..., overwrite = length(list(...)) > 0)
    return(colorVec)
}


#' @title Finish a SummarizedExperiment creation.
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment object
#' @importFrom viridis viridis
#' @importFrom SummarizedExperiment colData
#' @importFrom lubridate format_ISO8601
#' @noRd
finishExperiment <- function(exp) {

    exp <- exp[, order(exp$injection_time, method = "radix")]
    exp$datetime <- lubridate::as_datetime(exp$injection_time)


    exp$color <- setSampleColors(exp$type)[exp$type]

    if (!"batch" %in% colnames(colData(exp))) {
        exp$batch <- 1
    } else {
        exp$batch <- as.integer(as.factor(exp$batch))
    }

    if (!"injection" %in% colnames(colData(exp))) {
      colData(exp)$injection <- 1
  }

  if (!"order" %in% colnames(colData(exp))) {

      order <- colData(exp) %>%
          as.data.frame() %>%
          arrange(.data$injection_time) %>%
          group_by(.data$batch) %>%
          mutate(order = 1:n()) %>%
          pull(order)

      colData(exp)$order <- order
  }


  exp <- calculateRatio(exp) %>%
    identifyOutliers() %>%
    identifyMisInjections() %>%
    #addLinearRange(calType = "CAL", saveAssay = "CalRange") %>%
    addLinearRange(calType = "ACAL", saveAssay = "ACALRange")

  return(exp)
}

#' @title Fill in missing aliquot-compound combinations
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param df Dataframe in long format
#' @param rowIndex Column name in `df` that should be used for compounds.
#' @param colIndex Column name in `df` that should be used for aliquots.
#' @param rowData Dataframe with rowData
#' @param colData Dataframe with colData
#' @param assayName Column name in `df` that represents the area / intensities
#' of measured compounds.
#' @importFrom dplyr bind_rows
fillMissingCombinations <- function(df, rowIndex, colIndex, rowData, colData,
                                    assayNames) {

  comps <- names(which(table(df[, rowIndex]) < nrow(colData)))
  aliqs <- names(which(table(df[, colIndex]) < nrow(rowData)))

  a <- paste0(df[, colIndex], df[, rowIndex])

  grid <- expand.grid(a = aliqs, b = comps)
  b <- paste0(grid$a, grid$b)

  grid <- grid[!b %in% a, ]
  colnames(grid) <- c(colIndex, rowIndex)

  for (name in assayNames) {
    grid[[name]] <- NA
  }

  a <- dplyr::bind_rows(df[, c(colIndex, rowIndex, assayNames)], grid)
  a[order(a[, colIndex], a[, rowIndex], method = "radix"), ]
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
#' @returns a data.frame in a generic, human-readable format
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
#' # Convert to a table format
#' expToCombined(exp, rowIndex = "Compound", colIndex = "Aliquot")
expToCombined <- function(exp, rowIndex = "Compound", colIndex = "Aliquot") {
  if (!validateExperiment(exp)) {
    stop("Invalid SummarizedExperiment")
  }

  mat <- as.matrix(assay(exp))
  x <- reshape2::melt(mat)[, -3]
  colnames(x) <- c(rowIndex, colIndex)

  assays <- do.call(cbind, lapply(assayNames(exp), function(assay){
    mat <- as.matrix(assay(exp, assay))
    reshape2::melt(mat)[, 3, drop = FALSE]
  }))
  colnames(assays) <- assayNames(exp)
  cols <- c(colnames(x),
            colnames(assays),
            colnames(colData(exp)),
            colnames(rowData(exp)))

  df <- as.data.frame(cbind(x, assays,
                            colData(exp)[x[, colIndex], ],
                            rowData(exp)[x[, rowIndex], ] ))
  df <- df[order(df[, rowIndex], df[, colIndex]), ]
  rownames(df) <- seq_len(nrow(df))
  colnames(df) <- cols

  df[, !is.na(colnames(df))]
}

#' @title Add known concentrations of calibration lines
#' @description mzQuality supports the calculation of absolute concentrations
#' given concentrations of Calibration lines (CAL) are present. Using this
#' function, a new assay in the experiment is created called 'Concentration',
#' which is filled with NA and the supplied concentrations.
#'
#' The concentrations should be supplied as a data frame with rows as the
#' names of Calibration lines and columns as compound names.
#' Optionally, the `filterComps` argument can be set to TRUE, which only keeps
#' compounds in the experiment with a known concentration in the CAL lines.
#' @param exp SummarizedExperiment to add the data to
#' @param df Dataframe with concentrations of CAL lines
#' @param filterComps Boolean value, should compounds without concentration
#' be removed?
#' @returns SummarizedExperiment with an added assay called "Concentrations".
#' The values of the assay are filled with the values provided in `df`, other
#' values will be `NA`.
#' @export
addConcentrations <- function(exp, df, filterComps = FALSE) {

    if (!validateExperiment(exp)) stop("Invalid Experiment")
    rownames(df) <- make.names(rownames(df))
    if (!validateCalibration(exp, df)) stop("Invalid concentration dataframe")

    comps <- intersect(rownames(exp), rownames(df))
    comps

    if (!all(colnames(exp) %in% colnames(df))) {
        find <- gregexec("([A-Za-z]+)([0-9]+)", colnames(df))
        matches <- regmatches(colnames(df), find)
        matches <- do.call(rbind, lapply(matches, t))
        idx1 <- exp$calno %in% as.integer(matches[, 3])
        idx2 <- exp$type %in% matches[, 2]
        aliqs <- colnames(exp)[idx1 & idx2]
        if (!any(aliqs %in% colnames(exp))) stop("No matching aliquots found")
    } else {
        aliqs <- colnames(exp)[colnames(exp) %in% colnames(df)]
    }



    if (length(unique(exp[, aliqs]$type)) > 1) {
        stop("More than one Calibration type is not supported at this time")
    }
    metadata(exp)$concentration <- unique(exp[, aliqs]$type)


    assay(exp, "concentration") <- matrix(
        NA, nrow = nrow(exp), ncol = ncol(exp), dimnames = dimnames(exp)
    )

    # Subset concentration when not all expected callines are present
    cols <- aliqs[aliqs %in% colnames(exp[, exp$batch == unique(exp$batch)[1]])]
    df <- df[, exp[, cols]$calno + 1]

    res <- do.call(cbind, lapply(unique(exp$batch), function(batch){
        cols <- intersect(aliqs, colnames(exp[, exp$batch == batch]))
        colnames(df) <- cols
        df
    }))



    if (length(comps) > 0 & length(aliqs) > 0) {
        res <- do.call(cbind, lapply(unique(exp$batch), function(batch){
            cols <- intersect(aliqs, colnames(exp[, exp$batch == batch]))
            colnames(df) <- cols
            df
        }))

        assay(exp, "concentration")[comps, aliqs] <- as.matrix(res[comps, aliqs])
    }



    columns <- exp$type == metadata(exp)$concentration & exp$calno == 0
    concs <- assay(exp[, columns], "concentration")
    if (ncol(concs) > 0) {
        concs[is.na(concs)] <- 0
        assay(exp[, columns], "concentration") <- concs
    }


    if (filterComps) {
        df <- df[rowSums(df, na.rm = TRUE) > 0, ]
        exp <- exp[intersect(rownames(exp), rownames(df)), ]
    }

    rowData(exp)$hasKnownConcentrations <- FALSE
    rowData(exp)$hasCalculatedConcentrations <- FALSE
    rowData(exp[comps, ])$hasKnownConcentrations <- TRUE

    message("Added concentrations")
    exp
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
#' # Calculate the ratio
#' exp <- calculateRatio(exp)
#'
#' # Show results
#' assay(exp, "ratio")
calculateRatio <- function(exp, assay1 = metadata(exp)$primary,
                           assay2 = metadata(exp)$secondary) {


  if (!metadata(exp)$hasIS) {
    df <- assay(exp, assay1)
  } else {
    df <- assay(exp, assay1) / assay(exp, assay2)
  }

    df[df <= 0] <- NA
    assay(exp, "ratio") <- df
  #
  # assay(exp, "ratio") <- t(apply(df, 1, function(x) {
  #   ifelse(x <= 0, min(x[x > 0], na.rm = TRUE), x)
  # }))

  exp
}
