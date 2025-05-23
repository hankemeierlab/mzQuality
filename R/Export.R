#' @title Create a dynamic mzQuality summary report
#' @description Generates an interactive HTML report summarizing the full
#'     experiment across all compounds. Includes plots such as heatmaps, PCA
#'     plots, and quality control (QC) visualizations to help evaluate batch
#'     effects, sample distributions, and analytical variation.
#' @details This report helps to assess the overall quality of the experiment
#'     and offers a visual overview of sample-level patterns across selected
#'     assays. The plots are rendered interactively (via Plotly) for easier
#'     exploration and interpretation.
#' @param exp SummarizedExperiment object containing the experiment data.
#' @param folder Character string specifying the output directory for the
#'     report.
#' @param output File name for the HTML output. Defaults to
#'     "mzquality-report.html".
#' @param plots Character vector of plots to include in the report. Supported:
#'     "Aliquot", "PCA", "QC".
#' @param assays Character vector specifying which assays in `exp` to include.
#'     Only existing assay names will be used.
#' @returns An interactive HTML report is saved to the specified `folder`.
#' @noRd
summaryReport <- function(
        exp, folder, output = "mzquality-report.html",
        plots = c("Aliquot", "PCA", "QC"),
        assays = c("ratio", "ratio_corrected")
) {
    stopifnot(isValidExperiment(exp))
    assays <- intersect(assays, assayNames(exp))

    .createReport(
        path = file.path(folder, "Plots"),
        parameters = list(exp = exp, plots = plots, assays = assays),
        template = "summaryTemplate.Rmd",
        output = output
    )
}

#' @title Create dynamic compound report(s)
#' @description Generates interactive HTML reports for individual compounds
#'     from a SummarizedExperiment object. Each report includes compound-
#'     specific plots across selected assays for quality control and visual
#'     inspection.
#' @details The function parallelizes report generation over the number of
#'     compounds using `multiApply()`, supporting multicore environments. Each
#'     compound report is stored as a separate HTML file named after the
#'     compound.
#' @param exp SummarizedExperiment object containing the experiment data.
#' @param folder Output folder to store the compound reports.
#' @param assays Character vector of assay names to visualize in the reports.
#'     Must exist in `exp`.
#' @param project Character string denoting the project name. Defaults to
#'     "mzQuality".
#' @returns One HTML report per compound is saved in the specified folder.
#' @noRd
compoundReports <- function(
        exp, folder, assays = c("ratio", "ratio_corrected"),
        project = "mzQuality"
) {
    stopifnot(isValidExperiment(exp))
    outputNames <- sprintf("%s.html", make.names(rownames(exp)))

    assays <- intersect(assays, assayNames(exp))
    for (i in seq_len(nrow(exp))) {
        .createReport(
            path = file.path(folder, "Plots", "Compounds"),
            parameters = list(
                exp = exp[i, ],
                assays = assays,
                plots = "",
                project = project,
                vendor = "Unknown"
            ),
            template = "compoundTemplate.Rmd",
            output = outputNames[i]
        )
    }
}

#' @title Create an Rmarkdown report
#' @description Renders an Rmarkdown template into a self-contained HTML
#'     report.
#' @details Used internally to render experiment- or compound-specific
#'     reports. Handles creation of temporary intermediate directories and
#'     cleanup after rendering.
#' @param path Output directory where the final HTML report is stored.
#' @param parameters Named list of parameters to pass to the Rmarkdown
#'     template.
#' @param template Name of the Rmarkdown template file (located in the
#' "template" directory of the package).
#' @param output Output file name (HTML).
#' @returns HTML file is rendered and saved to disk. Function returns NULL
#'     invisibly.
#' @importFrom rmarkdown render
#' @importFrom utils assignInNamespace
#' @noRd
.createReport <- function(path, parameters, template, output) {
    file <- system.file("templates", template, package = "mzQuality")

    render(
        input = file,
        params = parameters,
        output_file = output,
        quiet = TRUE,
        output_dir = path,
        output_format = "html_document",
        output_options = list(self_contained = TRUE),
        clean = TRUE
    )
}

#' @title Export rowData, colData, and Assays
#' @description Produces human-readable tables of rowData, colData, and assays
#'   from a SummarizedExperiment object. The rowData is stored as
#'   'Compounds.tsv', colData as 'Aliquots.tsv', and assays are stored under
#'   their respective assay names.
#' @details This function creates a folder named "Export" in the specified
#'   directory and writes the rowData, colData, and assays of the
#'   SummarizedExperiment object into separate tab-separated files.
#' @returns Creates files of rowData, colData, and assays stored in `exp`.
#' @param exp A SummarizedExperiment object.
#' @param folder A folder where the tables will be stored.
#' @importFrom utils write.table
#' @export
#' @examples
#' # Read example dataset
#' exp <- readRDS(system.file("extdata/data.RDS", package = "mzQuality"))
#'
#' # Perform analysis
#' exp <- doAnalysis(exp)
#'
#' # Export tables in the experiment
#' exportTables(exp[1:10, 1:10], folder = tempdir())
exportTables <- function(exp, folder) {

    folder <- file.path(folder, "Exports")
    if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    }

    write.table(
        x = expToCombined(exp),
        file = file.path(folder, "mzQualityInput.tsv"),
        sep = "\t",
        row.names = TRUE,
        quote = FALSE
    )

    x <- cleanDataframe(colData(exp))
    write.table(
        file = file.path(folder, "Aliquots.tsv"),
        row.names = TRUE,
        x = x, sep = "\t"
    )


    x <- cleanDataframe(rowData(exp))

    write.table(
        file = file.path(folder, "Compounds.tsv"),
        row.names = TRUE,
        x = x, sep = "\t"
    )

    for (assay in assayNames(exp)) {
        write.table(
            file = file.path(folder, sprintf("%s.tsv", assay)),
            row.names = TRUE, x = assay(exp, assay), sep = "\t"
        )
    }
}

#' @title Generate a summary of aliquot-level metadata
#' @description Creates a summary data frame of aliquot-level metadata from
#'     the provided experiment object.
#' @details This function extracts aliquot-level metadata from the
#'     `colData` of the experiment object, excluding columns specified in
#'     `pkg.env$colDataExclude`. It adds a column indicating whether each
#'     aliquot is selected based on the presence of outliers and rounds
#'     numeric columns to three decimal places.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @return A data frame summarizing aliquot-level metadata, including a
#'     column for selection status and rounded numeric values.
#' @noRd
.getAliquotSummary <- function(exp){
    df <- colData(exp)
    df <- df[, colnames(df) != "color", drop = FALSE]

    aliquotSummary <- cleanDataframe(df)

    if (!"outlier" %in% colnames(aliquotSummary)) {
        aliquotSummary$outlier <- FALSE
    }
    aliquotSummary$selected <- ifelse(aliquotSummary$outlier, "No", "Yes")
    aliquotSummary <- aliquotSummary[, colnames(aliquotSummary) != "use"]

    aliquotSummary <- aliquotSummary %>%
        mutate_if(is.numeric, round, digits = 3)

    return(aliquotSummary)
}

#' @title Generate a summary of compound-level metadata
#' @description Creates a summary data frame of compound-level metadata from
#'     the provided experiment object.
#' @details This function extracts compound-level metadata from the
#'     `rowData` of the experiment object, excluding columns specified in
#'     `pkg.env$rowDataExclude`. It calculates a confidence level for each
#'     compound based on background signal and RSD thresholds, adds a column
#'     indicating whether each compound is selected, and rounds numeric
#'     columns to three decimal places.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param backgroundPercent Numeric, the threshold for background signal
#'     percentage. Defaults to 40.
#' @param cautionRSD Numeric, the RSD threshold for cautionary reporting.
#'     Defaults to 15.
#' @param nonReportableRSD Numeric, the RSD threshold for non-reportable
#'     data. Defaults to 30.
#' @return A data frame summarizing compound-level metadata, including
#'     confidence levels, selection status, and rounded numeric values.
#' @noRd
.getCompoundSummary <- function(
        exp, backgroundPercent = 40, cautionRSD = 15, nonReportableRSD = 30
) {
    df <- rowData(exp)
    df <- df[!colnames(df) %in% pkg.env$rowDataExclude]

    compoundSummary <- cleanDataframe(df)

    # Use case_when for cleaner conditional logic
    compoundSummary$confidence <- dplyr::case_when(
        compoundSummary$rsdqcCorrected > nonReportableRSD ~ "Non Reportable",
        compoundSummary$backgroundSignal > backgroundPercent ~ "Low SNR",
        compoundSummary$backgroundSignal <= backgroundPercent &
            compoundSummary$rsdqcCorrected < cautionRSD ~ "High",
        TRUE ~ "Caution"
    )

    compoundSummary <- compoundSummary %>%
        mutate(selected = ifelse(.data$use, "Yes", "No")) %>%
        select(-use) %>%
        mutate_if(is.numeric, round, digits = 3)

    return(compoundSummary)
}

#' @title Write the new export file
#' @description Exports data from a SummarizedExperiment object to an Excel
#'     file, including summaries and assay data.
#' @details This function generates an Excel file containing multiple sheets
#'     with summaries and assay data from the provided SummarizedExperiment
#'     object. It includes:
#'     - A summary of aliquots, compounds, and internal standards.
#'     - Aliquot-level and compound-level metadata summaries.
#'     - Assay data, excluding those specified in `pkg.env$assayExclude`.
#'     - Additional metadata such as batch R-squared values, studentized
#'       residuals, and linear range fractions, if available.
#'
#'     The function handles Excel sheet name length limitations by truncating
#'     names to 31 characters.
#' @param folder A string specifying the directory where the Excel file will
#'     be saved.
#' @param file A string specifying the name of the Excel file to be created.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param types A character vector specifying the sample types to include in
#'     the export. Defaults to all types in `exp$type`.
#' @param backgroundPercent Numeric, the threshold for background signal
#'     percentage. Defaults to 40.
#' @param cautionRSD Numeric, the RSD threshold for cautionary reporting.
#'     Defaults to 15.
#' @param nonReportableRSD Numeric, the RSD threshold for non-reportable data.
#'     Defaults to 30.
#' @param selectedOnly Logical, whether to include only selected aliquots and
#'     compounds. Defaults to `FALSE`.
#' @return None. The function writes an Excel file as a side effect.
#' @importFrom dplyr mutate_if
#' @importFrom openxlsx write.xlsx
#' @importFrom SummarizedExperiment assays
#' @noRd
.writeExcelExport <- function(
        folder, file, exp, types = exp$type,
        backgroundPercent = 40, cautionRSD = 15,
        nonReportableRSD = 30, selectedOnly = FALSE
) {

    IS <- unique(rowData(exp)$compound_is)
    overallSummary <- data.frame(
        Info = c(ncol(exp), nrow(exp), length(IS)),
        row.names = c("Aliquots", "Compounds", "Internal Standards")
    )

    if (selectedOnly) {
        exp <- exp[rowData(exp)$use, exp$use]
    }

    exp <- exp[, exp$type %in% types]

    tabs <- list(
        Summary = overallSummary,
        Aliquot = .getAliquotSummary(exp),
        Compound = .getCompoundSummary(exp)
    )


    assayList <- assays(exp)[!assayNames(exp) %in% pkg.env$assayExclude]
    isNumeric <- vapply(assayList, function(m) all(is.numeric(m)), logical(1))
    assayList <- assayList[isNumeric]

    tabs <- c(tabs, assayList[sort(names(assayList))])

    if (!is.null(metadata(exp)$concentration)) {
        tabs <- c(tabs, list(
            BatchR2 = rowData(exp)$concentrationR2,
            StudentizedResiduals = rowData(exp)$studentizedResiduals,
            LinearRangeFraction = rowData(exp)$linearRanges
        ))
    }

    # Excel limitations
    names(tabs) <- substr(names(tabs), 0, 31)

    path <- file.path(folder, "Exports", file, fsep = "/")
    write.xlsx(
        x = tabs, file = path,
        rowNames = TRUE, colNames = TRUE,
        keepNA = TRUE, na.string = "NA"
    )
}

#' @title Create a project-specific report folder structure
#' @description Creates a folder structure for storing reports, plots, and
#'     compound-specific data for a given project.
#' @details This function generates a folder structure within a specified base
#'     folder for a given project. The structure includes:
#'     - A main project folder
#'     - A subfolder for exports
#'     - A subfolder for plots
#'     - A subfolder for compound-specific plots
#'
#'     If the folders do not already exist, they are created recursively.
#' @param baseFolder A string specifying the base directory where the project
#'     folder will be created.
#' @param projectName A string specifying the name of the project. This will
#'     be used as the name of the main project folder.
#' @return None. The function creates the folder structure as a side effect.
#' @noRd
.createReportFolder <- function(outputFolder){

    if (!dir.exists(outputFolder)) {
        dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
    }

    reports <- file.path(outputFolder, "Exports")
    plotFolder <- file.path(outputFolder, "Plots")
    compoundFolder <- file.path(plotFolder, "Compounds")

    for (folder in c(reports, plotFolder, compoundFolder)) {
        if (!dir.exists(folder)) {
            dir.create(folder, recursive = TRUE, showWarnings = FALSE)
        }
    }
}

#' @title Generate project reports and export data
#' @description Creates a folder structure, exports data tables, and generates
#'     summary and compound-specific reports for a given project.
#' @details This function automates the process of creating reports for a
#'     project. It performs the following tasks:
#'     - Creates a folder structure for storing reports and plots.
#'     - Exports data tables to the specified folder.
#'     - Generates a summary report with specified plots.
#'     - Generates compound-specific reports.
#'
#'     The function allows customization of the report content and thresholds
#'     for background percentage and relative standard deviation (RSD).
#' @param folder A string specifying the base directory where the reports will
#'     be saved.
#' @param project A string specifying the name of the project. This will be
#'     used to create a project-specific folder.
#' @param exp A SummarizedExperiment object containing the experimental data.
#' @param summaryPlots A character vector specifying the types of summary plots
#'     to include in the summary report. Defaults to
#'     `c("Aliquot", "PCA", "QC")`.
#' @param makeSummaryReport Logical, whether to generate a summary report.
#'     Defaults to `TRUE`.
#' @param makeCompoundReport Logical, whether to generate compound-specific
#'     reports. Defaults to `TRUE`.
#' @param backgroundPercent Numeric, the percentage threshold for background
#'     correction. Defaults to `40`.
#' @param cautionRSD Numeric, the RSD threshold for cautionary reporting.
#'     Defaults to `15`.
#' @param nonReportableRSD Numeric, the RSD threshold for non-reportable data.
#'     Defaults to `30`.
#' @param assays A character vector specifying the assay names to include in
#'     the reports. Defaults to `c("ratio", "ratio_corrected")`.
#' @return None. The function creates reports and exports data as a side effect.
#' @examples
#' exp <- readRDS(system.file("extdata/data.RDS", package = "mzQuality"))
#'
#' exp <- doAnalysis(exp, removeOutliers = TRUE)
#'
#' createReports(
#'     folder = tempdir(),
#'     project = "MyProject",
#'     exp = exp[1:2, ],
#'     summaryPlots = c("Aliquot", "PCA"),
#'     makeSummaryReport = TRUE,
#'     makeCompoundReport = TRUE,
#'     backgroundPercent = 50,
#'     cautionRSD = 10,
#'     nonReportableRSD = 25,
#'     assays = "ratio"
#' )
#' @export
createReports <- function(
        folder, project, exp, summaryPlots = c("Aliquot", "PCA", "QC"),
        makeSummaryReport = TRUE, makeCompoundReport = TRUE,
        backgroundPercent = 40, cautionRSD = 15, nonReportableRSD = 30,
        assays = c("ratio", "ratio_corrected")
) {

    folder <- file.path(folder, project, fsep = "/")
    # Create Report folder if it doesn't exist
    .createReportFolder(folder)
    exportTables(exp, folder)

    .writeExcelExport(
        folder = folder, file = "FinalReport.xlsx", exp = exp,
        backgroundPercent = backgroundPercent, cautionRSD = cautionRSD,
        nonReportableRSD = nonReportableRSD
    )

    if (makeSummaryReport) {
        summaryReport(
            exp = exp, folder = folder,
            plots = summaryPlots, assays = assays
        )
    }

    if (makeCompoundReport) {
        compoundReports(exp = exp, folder = folder, assays = assays)
    }
}
