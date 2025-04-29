library(testthat)
suppressPackageStartupMessages(library(mzQuality))
# Read the example dataset
dataFile <- system.file("data.RDS", package = "mzQuality")
exp <- readRDS(dataFile)


test_that("Internal Standards can be replaced", {
    exp2 <- replaceInternalStandards(exp)
    expect_true(is(exp, "SummarizedExperiment"))

    expect_true(all(rowData(exp2)$suggestedIS == rowData(exp2)$compound_is))
})
