#' @title Export Excel workbook
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment object
#' @param outfile File to export to
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom stats setNames
#' @returns A workbook object from `openxlsx`
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
#' # Export used aliquots and metabolites
#' exportFilterWorkbook(exp, outfile = tempfile())
exportFilterWorkbook <- function(exp, outfile = "Metabolites_to_keep.xlsx") {

    metabolites <- data.frame(
        Compound = rownames(exp),
        Used = rowData(exp)$use
    )

    aliquots <- data.frame(
        Aliquot = colnames(exp),
        Used = exp$use
    )

    tabs <- list(
        Summary = summaryTable(exp),
        Metabolites = metabolites,
        Aliquots = aliquots
    )

    openxlsx::write.xlsx(
        x = tabs,
        file = outfile,
        rowNames = TRUE,
        colNames = TRUE,
        keepNA = TRUE,
        na.string = "NA"
    )
}


#' @title Format SummarizedExperiment to dataframe
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment object
#' @noRd
formatExport <- function(exp) {
    exp <- exp[order(rownames(exp)), order(colnames(exp))]
    ratios <- t(as.data.frame(assay(exp, "ratio_corrected")))
    rsdqc <- t(data.frame(RSDQC = rowData(exp)$rsdqcCorrected))

    if ("backgroundSignal" %in% colnames(rowData(exp))) {
        be <- t(data.frame(Background = rowData(exp)$backgroundSignal * 100))
    } else {
        be <- t(data.frame(Background = rep(NA, nrow(exp))))
    }

    be[is.nan(be)] <- NA
    rsdqc[is.nan(rsdqc)] <- NA
    ratios[is.nan(ratios)] <- NA

    colnames(be) <- rownames(exp)
    colnames(rsdqc) <- rownames(exp)
    rsdqc <- cbind(Type = NA, Batch = NA, rbind(rsdqc, be))
    as.data.frame(rbind(rsdqc, cbind(Type = exp$type, Batch = exp$batch, ratios)))
}



#' @title Retrieve a summary of the SummarizedExperiment
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment object
#' @param project Project name, defaults to "mzQuality"
#' @param vendor Vendor name of the files, defaults to "Unknown"
#' @returns data.frame with a column "Info" and several rows of basic
#' information.
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
#' # Request a summary in table format
#' summaryTable(exp, project = "MyProject", vendor = "Sciex")
summaryTable <- function(exp, project = "mzQuality", vendor = "Unknown") {
    df <- t(data.frame(
        Project = project,
        Date = lubridate::format_ISO8601(metadata(exp)$Date),
        Vendor = vendor,
        QC = metadata(exp)$QC,
        Batches = length(unique(exp$batch)),
        Aliquots = ncol(exp), Compounds = nrow(exp),
        Version = as.character(utils::packageVersion("mzQuality2"))
    ))
    colnames(df) <- "Info"
    as.data.frame(df)
}

#' @title Export to Excel
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment Object
#' @param file Name of the file to export to
#' @importFrom openxlsx write.xlsx
#' @returns A workbook object from `openxlsx`
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
#' # Request a summary in table format
#' exportExcel(exp, file = tempfile())
exportExcel <- function(exp, file,  backgroundPercent = 40, cautionRSD = 15, nonReportableRSD = 30) {
    exp <- exp[, exp$use]

    df <- summaryTable(exp)

    qcPresenceThreshold <- rowData(exp)$qcPresenceThreshold

    hc <- filterRSDQC(exp[qcPresenceThreshold, ], min = 0, max = cautionRSD)
    if (all(dim(hc) > 0)) {
        hc <- filterBackground(hc, include.na = TRUE, max = backgroundPercent / 100)
    }
    hc <- formatExport(hc)

    lc <- filterRSDQC(exp[qcPresenceThreshold, ], min = cautionRSD, max = nonReportableRSD)
    if (all(dim(lc) > 0)) {
        lc <- filterBackground(lc, include.na = TRUE, max = backgroundPercent / 100)
    }
    lc <- formatExport(lc)
    comps <- exp[setdiff(rownames(exp), c(colnames(hc), colnames(lc))), ]



    non_reported <- formatExport(comps)


    openxlsx::write.xlsx(list(
        Summary = df, `High Confidence` = hc,
        `With Caution` = lc, `Low SignalNoise` = non_reported
    ),
    file = file,
    rowNames = TRUE, colNames = TRUE, keepNA = TRUE, na.string = "NA"
    )
}


#' @title Create a dynamic mzQuality report
#' @description A summary report is a dynamic report containing plots involving
#' all compounds. Examples include a heatmap, aliquot plot and violin plot.
#' All plots are interactive and can be used for assessing experiments at a
#' later stage.
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment Object
#' @param folder Directory where to store the report.
#' @param project Project name, defaults to "mzQuality"
#' @param vendor Vendor name of the files, defaults to "Unknown"
#' @param output Output file
#' @returns Creates a HTML file with dynamic visualizations
#' @export
#' @examples
#' # Read example dataset
#' data(mzQualityExp)
#'
#' # Request a summary report
#' # summaryReport(mzQualityExp,
#' #    folder = tempdir(), project = "MyProject",
#' #    vendor = "Agilent", output = tempfile()
#' #)
summaryReport <- function(exp, folder, output = "mzquality-report.html",
                          plots = c("Aliquot", "PCA", "QC"),
                          assays = c("ratio", "ratio_corrected")) {

    if (!validateExperiment(exp)) stop("Invalid Experiment")

    assays <- assays[assays %in% assayNames(exp)]

    create_report(
        path = file.path(folder),
        parameters = list(exp = exp, plots = plots, assays = assays),
        template = "mzquality-report-se.Rmd",
        output = output
    )
}

#' @title Report job
#' @noRd
compoundReportJob <- function(i, exp, folder, assays, project, vendor) {
  output <- glue::glue("{comp}.html", comp = make.names(rownames(exp)[i]))
  create_report(
    path = file.path(folder),
    parameters = list(
      exp = exp[i, ],
      assays = assays,
      plots = "",
      project = project,
      vendor = vendor
    ),
    template = "mzquality-compounds-se.Rmd",
    output = output
  )
  NULL
}

#' @title Create dynamic compound report(s)
#' @description A compound report is a dynamic report containing plots
#' regarding a single compound. By default, this function produces a plot for
#' each compound, but can also be made for a single compound.
#' @details placeholder
#' @returns placeholder
#' @param exp SummarizedExperiment object
#' @param folder Folder where to store the compound report.
#' @param compounds Which compounds to create a report for. Defaults to all
#' compounds in `exp`
#' @param project Project name, defaults to "mzQuality"
#' @param vendor Vendor name of the files, defaults to "Unknown"
#' @param output Output file
#' @returns Creates a HTML file with dynamic visualizations
#' @importFrom foreach foreach %dopar% %do%
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
#' # Request a compound report
compoundReports <- function(exp, folder, assays = c("ratio", "ratio_corrected"),
                            project = "mzQuality", cores = 1, verbose = F, shinyProgress = NULL) {

    if (!validateExperiment(exp)) stop("Invalid Experiment")

    multiApply(
      func = compoundReportJob,
      jobs = nrow(exp),
      cores = cores,
      shinyProgress = shinyProgress,
      progressTitle = "Construction Compound reports",
      verbose = verbose,
      exp = exp,
      folder = folder,
      assays = assays,
      project = project,
      vendor = "Unknown"
    )
}

#' @title Create an Rmarkdown report
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param path Where to store the report
#' @param parameters List of parameters
#' @param template Name of template file
#' @param output Name of the report
#' @noRd
create_report <- function(path, parameters, template, output) {
    file <- system.file("rmd", template, package = "mzQuality2")

    library(rmarkdown)

    assignInNamespace("clean_tmpfiles", function() {}, ns = "rmarkdown")

    intermed <- paste0("intermediates_", stringi::stri_rand_strings(1, 10))
    rmarkdown::render(
      input = file,
      params = parameters,
      output_file = output,
      quiet = FALSE,
      output_dir = path,
      output_format = "html_document",
      output_options = list(self_contained = T),
      intermediates_dir = intermed,
      clean = T
    )
    system(glue::glue("rm -r {intermed}"))
}

#' @title Export rowData, colData, and assays
#' @description This function produces human-readable tables of rowData,
#' colData and assays of a SummarizedExperiment object. rowData will be stored
#' as 'Compounds.tsv', colData will be stored as 'Aliquots.tsv' and assays will
#' be stored under their assay name.
#' @details placeholder
#' @returns placeholder
#' @param exp A SummarizedExperiment object
#' @param folder Folder where to store the tables
#' @returns Creates files of rowData, colData and assays stored in `exp`
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
#' # Do Analysis
#' exp <- doAnalysis(exp)
#'
#' # Export tables in Experiment
#' exportTables(exp, folder = tempdir())
exportTables <- function(exp, folder) {
    if (!validateExperiment(exp)) stop("Invalid Experiment")

    folder <- file.path(folder)


    cols <- sapply(colData(exp), function(x) !"matrix" %in% class(x))

    x <- colData(exp)[, cols]
    utils::write.table(
        file = file.path(folder, "Aliquots.tsv"),
        row.names = TRUE,
        x = as.data.frame(x), sep = "\t"
    )


    cols <- sapply(rowData(exp), function(x) !"matrix" %in% class(x))
    x <- rowData(exp)[, cols]

    utils::write.table(
        file = file.path(folder, "Compounds.tsv"),
        row.names = TRUE,
        x = as.data.frame(x), sep = "\t"
    )

    for (assay in assayNames(exp)) {
        utils::write.table(
            file = file.path(folder, sprintf("%s.tsv", assay)),
            row.names = TRUE, x = assay(exp, assay), sep = "\t"
        )
    }
}


#' @title Zip the folder with results
#' @description placeholder
#' @details placeholder
#' @returns placeholder
#' @param file output name of file
#' @param output Folder to be zipped
#' @returns A zip file is created of the folders created with the export
#' functions and reports.
#' @export
zipFolder <- function(file, output) {
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zip(zipfile = file, files = list.files(output), root = output)
  }
}

#' @title Write the new export file
#' @description
#' @details
#' @returns
#' @export
writeNewExport <- function(outFile, exp, types = exp$type,
                           backgroundPercent = 40, cautionRSD = 15, nonReportableRSD = 30,
                           digits = 3, selectedOnly = FALSE){



    overallSummary <- data.frame(
        Info = c(
            ncol(exp), nrow(exp), length(unique(rowData(exp)$compound_is))
        )
    )
    rownames(overallSummary) <- c(
        "Aliquots",
        "Compounds",
        "Internal Standards"
    )

    if (selectedOnly) {
        exp <- exp[rowData(exp)$use, exp$use]
    }

    toUse <- which(!colnames(colData(exp)) %in% pkg.env$colDataExclude)
    aliquotSummary <- as.data.frame(colData(exp)[, toUse])
    aliquotSummary$selected <- ifelse(aliquotSummary$outlier, "No", "Yes")
    aliquotSummary <- aliquotSummary[, colnames(aliquotSummary) != "use"]

    toUse <- which(!colnames(rowData(exp)) %in% pkg.env$rowDataExclude)
    compoundSummary <- as.data.frame(rowData(exp)[, toUse])

    compoundSummary$confidence <- ifelse(compoundSummary$backgroundSignal <= backgroundPercent & compoundSummary$rsdqcCorrected < cautionRSD, "High", "Caution")
    compoundSummary$confidence <- ifelse(compoundSummary$rsdqcCorrected > nonReportableRSD, "Non Reportable", compoundSummary$confidence)
    compoundSummary$confidence <- ifelse(compoundSummary$backgroundSignal > backgroundPercent, "Low SNR", compoundSummary$confidence)
    compoundSummary$selected <- ifelse(compoundSummary$use, "Yes", "No")
    compoundSummary <- compoundSummary[, colnames(compoundSummary) != "use"]

    compoundSummary <- compoundSummary %>%
        mutate_if(is.numeric, round, digits = digits)

    aliquotSummary <- aliquotSummary %>%
        mutate_if(is.numeric, round, digits = digits)

    exp <- exp[, exp$type %in% types]


    assays(exp) <- assays(exp)[!assayNames(exp) %in% pkg.env$assayExclude]
    assayList <- lapply( assays(exp), function(m){
        if (all(is.numeric(m))) {
            m <- round(m, digits)
        }
        return(cbind(Batch = exp$batch, t(m)))
    })

    names(assayList) <- assayNames(exp)

    tabs <- list(
        Summary = overallSummary,
        Aliquot = aliquotSummary,
        Compound = compoundSummary
    )
    tabs <- c(tabs, assayList[sort(names(assayList))])

    if (!is.null(metadata(exp)$concentration)) {
        tabs <- c(tabs, list(
            BatchR2 = rowData(exp)$concentrationR2,
            StudentizedResiduals = rowData(exp)$studentizedResiduals,
            LinearRangeFraction = rowData(exp)$linearRanges
        ))
    }

    names(tabs) <- substr(names(tabs), 0, 31)

    openxlsx::write.xlsx(
        x = tabs,
        file = outFile,
        rowNames = TRUE,
        colNames = TRUE,
        keepNA = TRUE,
        na.string = "NA"
    )
}

#' @title Create Concentration report
#' @description
#' @details
#' @returns
#' @export
reportConcentrationsPerBatch <- function(exp, file){
    highConfLinearRange <- rowData(exp)$linearRanges > 0.7
    highR2 <- rowData(exp)$concentrationR2 > 0.95

    codes <- matrix(nrow = nrow(highR2), ncol = ncol(highR2), dimnames = dimnames(highR2))
    codes[highConfLinearRange & highR2] <- 1
    codes[highConfLinearRange & !highR2] <- 2
    codes[!highConfLinearRange & highR2] <- 3
    codes[!highConfLinearRange & !highR2] <- 4




    sheets <- lapply(unique(exp$batch), function(batch){
        x <- exp[, exp$batch == batch]

        header <- t(data.frame(
            Sample.Name = rownames(x),
            RSDQC = rowData(x)$rsdqcCorrected,
            Background = rowData(x)$backgroundSignal,
            Concentration = codes[, batch]
        ))
        rbind(
            header,
            t(assay(x, "concentration"))
        )
    })

    names(sheets) <- glue::glue("Batch {batch}", batch = unique(exp$batch))

    openxlsx::write.xlsx(
        x = sheets,
        file = file,
        rowNames = TRUE,
        keepNA = TRUE,
        na.string = "NA"
    )
}

#' @title Download zip file
#' @description
#' @details
#' @returns
#' @param exp
#' @export
downloadZip <- function(project, exp,
                        summaryPlots = c("Aliquot", "PCA", "QC"),
                        summaryReport = TRUE, compoundReport = TRUE,
                        backgroundPercent = 40, cautionRSD = 15, nonReportableRSD = 30,
                        copyDataLake = FALSE,
                        assays = c("ratio", "ratio_corrected"), cores = 1, shinyProgress = NULL) {

    output_folder <- file.path(getwd(), project)

    reports <- file.path(output_folder, "Exports")
    plotFolder <- file.path(output_folder, "Plots")
    compoundFolder <- file.path(plotFolder, "Compounds")

    dir.create(reports, recursive = TRUE, showWarnings = FALSE)
    dir.create(compoundFolder, recursive = TRUE, showWarnings = FALSE)

    message("Exporting Tables")
    exportTables(exp, reports)
    message("Exported Tables")

    timeFormat <- format(lubridate::now(), "%d-%m-%Y %H.%M.%S")
    file <- glue::glue("mzQuality_Original_{time}.tsv", time = timeFormat)

    utils::write.table(
        x = expToCombined(exp),
        file =  file.path(reports, file),
        sep = "\t",
        row.names = TRUE,
        quote = FALSE
    )
    message("Exported Original File")

    file <- "FinalReport_All_v2.xlsx"
    writeNewExport(
        file.path(reports, file),
        exp = exp,
        backgroundPercent = backgroundPercent,
        cautionRSD = cautionRSD,
        nonReportableRSD = nonReportableRSD
    )

    file <- "FinalReport_Sample_SelectedComps_v2.xlsx"
    writeNewExport(
        file.path(reports, file),
        exp = exp[rowData(exp)$use, ],
        types = "SAMPLE",
        backgroundPercent = backgroundPercent,
        cautionRSD = cautionRSD,
        nonReportableRSD = nonReportableRSD
    )

    # file <- "NonReported_QC_corrected_all_final_report.xlsx"
    # exportExcel(x, file.path(reports, file))

    message("Exporting Filter Workbook")
    file <- "Metabolites_to_keep.xlsx"
    exportFilterWorkbook(exp, file.path(reports, file))

    message("Exporting QC Corrected Workbook for Samples")

    file <- "QC_corrected_sample_final_report.xlsx"
    exportExcel(
        exp[, exp$type == "SAMPLE"],
        file.path(reports, file),
        backgroundPercent = backgroundPercent,
        cautionRSD = cautionRSD,
        nonReportableRSD = nonReportableRSD
    )

    message("Exporting QC Corrected Workbook for All")

    file <- "QC_corrected_all_final_report.xlsx"
    exportExcel(
        exp,
        file.path(reports, file),
        backgroundPercent = backgroundPercent,
        cautionRSD = cautionRSD,
        nonReportableRSD = nonReportableRSD
    )


    if (summaryReport) {
        message("Generating Summary Report")
        summaryReport(exp = exp, folder = plotFolder, plots = summaryPlots, assays = assays)
        message("Created Summary Report")
    }

    if (compoundReport) {
        message("Creating Compound Reports")
        compoundReports(exp, folder = compoundFolder, assays = assays, cores = cores, shinyProgress = shinyProgress)
        message("Created Compound Reports")
    }
}

