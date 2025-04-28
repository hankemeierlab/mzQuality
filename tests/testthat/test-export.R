library(testthat)
suppressPackageStartupMessages(library(mzQuality))
# Read the example dataset
combinedFile <- system.file("example.tsv", package = "mzQuality")

dataFile <- system.file("data.RDS", package = "mzQuality")
concentrationsFile <- system.file("concentrations.txt", package = "mzQuality")

exp <- readRDS(dataFile)
concentrations <- read.delim(concentrationsFile)

test_that("A combined file can be converted", {
    combined <- buildCombined(combinedFile)
    x <- buildExperiment(combined)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(validateExperiment(x))
})

test_that("exportFilterWorkbook works", {
    f <- tempfile()
    exportFilterWorkbook(exp, f)
    expect_true(file.exists(f))
})

test_that("formatExport works", {
    df <- mzQuality:::formatExport(exp)
    expect_true(is(df, "data.frame"))
    expect_true(all(c("Type", "Batch") %in% colnames(df)))
})

test_that("summaryTable  works", {
    df <- summaryTable(exp)
    expect_true(is(df, "data.frame"))
})

test_that("exportExcel  works", {
    f <- tempfile()
    exportExcel(exp, f)

    expect_true(file.exists(f))
})

test_that("summaryReport  works", {
    dir <- tempdir()

    summaryReport(exp[1:10, exp$type != "SAMPLE"], dir)

    path <- file.path(dir, "mzquality-report.html")
    expect_true(file.exists(path))
})

test_that("compoundReports   works", {
    dir <- tempdir()
    compoundReports(exp[1, ], dir, verbose = T)
    path <- file.path(dir, paste0(rownames(exp)[1], ".html"))
    expect_true(file.exists(path))
})

test_that("exportTables works", {
    dir <- tempdir()
    exportTables(exp, dir)

    expect_true(file.exists(file.path(dir, "Aliquots.tsv")))
    expect_true(file.exists(file.path(dir, "Compounds.tsv")))
})
