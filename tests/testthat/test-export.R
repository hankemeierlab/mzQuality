library(testthat)
suppressPackageStartupMessages(library(mzQuality))
# Read the example dataset
combinedFile <- system.file("example.tsv", package = "mzQuality")

dataFile <- system.file("data.RDS", package = "mzQuality")
concentrationsFile <- system.file("concentrations.txt", package = "mzQuality")

exp <- readRDS(dataFile)
concentrations <- read.delim(concentrationsFile)

test_that("A combined file can be converted", {
    combined <- readData(combinedFile)
    x <- buildExperiment(combined)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(validateExperiment(x))
})
