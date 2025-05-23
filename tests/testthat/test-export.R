library(testthat)
suppressPackageStartupMessages(library(mzQuality))
# Read the example dataset
combinedFile <- system.file("extdata/example.tsv", package = "mzQuality")

dataFile <- system.file("extdata/data.RDS", package = "mzQuality")
concentrationsFile <- system.file("extdata/concentrations.tsv", package = "mzQuality")

exp <- readRDS(dataFile)
concentrations <- read.delim(concentrationsFile)

exp <- doAnalysis(exp, removeOutliers = TRUE)

test_that("A combined file can be converted", {
    combined <- readData(combinedFile)
    x <- buildExperiment(combined)
    expect_true(isValidExperiment(x))
})


test_that("exportTables creates expected tsv files", {
    temp_folder <- tempdir()

    exportTables(exp[1:10, 1:10], folder = temp_folder)

    expect_true(file.exists(file.path(temp_folder, "Exports", "Aliquots.tsv")))
    expect_true(file.exists(file.path(temp_folder, "Exports", "Compounds.tsv")))
})


test_that("createReports generates reports and exports data", {
    temp_folder <- tempdir()

    subExp <- exp[1:2, exp$batch == exp$batch[1] & exp$type %in% c("SQC", "LQC")]

    createReports(
        folder = temp_folder,
        project = "TestProject",
        exp = subExp,
        summaryPlots = c("Aliquot", "PCA"),
        makeSummaryReport = TRUE,
        makeCompoundReport = TRUE,
        backgroundPercent = 50,
        cautionRSD = 10,
        nonReportableRSD = 25,
        assays = c("ratio")
    )

    expect_true(dir.exists(file.path(temp_folder, "TestProject")))
    expect_true(file.exists(file.path(temp_folder, "TestProject", "Exports", "FinalReport.xlsx")))
})

