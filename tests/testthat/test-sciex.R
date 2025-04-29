library(testthat)
suppressPackageStartupMessages(library(mzQuality))
# Read the example dataset
dataFile <- system.file("sciex.txt", package = "mzQuality")


test_that("Combined file can be constructed from sciex os file", {
    df <- readData(dataFile)
    expect_true(is(df, "data.frame"))
    expect_true(validateDataframe(df))
})

test_that("Combined file can be constructed from multiple sciex os files", {
    df <- readData(c(dataFile, dataFile, dataFile))
    expect_true(is(df, "data.frame"))
    expect_true(validateDataframe(df))
})
