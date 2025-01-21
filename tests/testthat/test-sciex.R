library(testthat)
suppressPackageStartupMessages(library(mzQuality2))
# Read the example dataset
dataFile <- system.file("sciex.txt", package = "mzQuality2")


test_that("Combined file can be constructed from sciex os file", {
    df <- buildCombined(dataFile)
    expect_true(is(df, "data.frame"))
    expect_true(validateDataframe(df))
})

test_that("Combined file can be constructed from multiple sciex os files", {
    df <- buildCombined(c(dataFile, dataFile, dataFile))
    expect_true(is(df, "data.frame"))
    expect_true(validateDataframe(df))
})
