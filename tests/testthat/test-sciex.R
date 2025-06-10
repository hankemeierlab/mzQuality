library(testthat)
suppressPackageStartupMessages(library(mzQuality))
path <- system.file("extdata/example.tsv", package = "mzQuality")
dataFile <- system.file("extdata/sciex.tsv", package = "mzQuality")


test_that("Combined file can be constructed from sciex os file", {
    df <- readData(dataFile)
    expect_true(is(df, "data.frame"))
    expect_true(isValidDataframe(df))
})

test_that("Combined file can be constructed from multiple sciex os files", {
    df <- readData(c(dataFile, dataFile, dataFile))
    expect_true(is(df, "data.frame"))
    expect_true(isValidDataframe(df))
})
