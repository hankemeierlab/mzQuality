library(testthat)
suppressPackageStartupMessages(library(mzQuality2))
# Read the example dataset
exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))


test_that("mzQuality export can be made", {
    file <- tempfile()
    exportFilterWorkbook(exp, outfile = file)
    expect_true(file.exists(file))

    file.remove(file)

    exportExcel(exp, file)
    expect_true(file.exists(file))
    file.remove(file)
})

test_that("New Export file can be made", {
    file <- tempfile()
    exp <- doAnalysis(exp, doAll = TRUE)
    writeNewExport(file, exp, types = exp$type,
                   backgroundPercent = 40, cautionRSD = 15, nonReportableRSD = 30,
                   digits = 3, selectedOnly = FALSE)
    expect_true(file.exists(file))
})

test_that("Concentration per batchcan be exported", {
    file <- tempfile()
    reportConcentrationsPerBatch(exp, file)
    expect_true(file.exists(file))
})
