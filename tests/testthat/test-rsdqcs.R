library(testthat)
suppressPackageStartupMessages(library(mzQuality))
path <- system.file("extdata", "example.tsv", package = "mzQuality")
exp <- buildExperiment(readData(path))

test_that("Internal Standards can be replaced", {
    exp2 <- replaceInternalStandards(exp)
    expect_true(isValidExperiment(exp))
    expect_true(all(rowData(exp2)$suggestedIS == rowData(exp2)$compound_is))
})
