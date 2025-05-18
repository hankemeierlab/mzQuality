library(testthat)
suppressPackageStartupMessages(library(mzQuality))
# Read the example dataset
dataFile <- system.file("data.RDS", package = "mzQuality")
concentrationsFile <- system.file("concentrations.txt", package = "mzQuality")

exp <- readRDS(dataFile)
concentrations <- read.delim(concentrationsFile)

file <- system.file("example.tsv", package = "mzQuality")
concentrationsFile <- system.file("concentrations.txt", package = "mzQuality")

exp <- doAnalysis(exp, doAll = TRUE, removeOutliers = TRUE)



test_that("A tab-delimited file can be converted into a SummarizedExperiment", {
    combined <- readData(file)
    expect_true(is(combined, "data.frame"))

    exp <- buildExperiment(combined)
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(isValidExperiment(exp))

    # Check that the last row is not included in the experiment
    # This should trigger fix missing
    exp <- buildExperiment(combined[-nrow(combined), ])
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(isValidExperiment(exp))
})


test_that("Is a valid experiment", {
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(isValidExperiment(exp))
    expect_true("area" %in% assayNames(exp))
    expect_true("ratio" %in% assayNames(exp))
    expect_true("ratio_corrected" %in% assayNames(exp))
})

test_that("experiment can be converted to combined", {
    comb <- expToCombined(exp)
    expect_true(is(comb, "data.frame"))
    expect_equal(nrow(comb), nrow(exp) * ncol(exp))
    expect_true(all(c("aliquot", "compound", "area", "injection_time", "type", "batch") %in% colnames(comb)))
})


test_that("Concentrations can be added", {
    x <- addConcentrations(exp, concentrations)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
    expect_true("hasKnownConcentrations" %in% colnames(rowData(exp)))

    expect_true(sum(rowData(exp)$hasKnownConcentrations) == nrow(concentrations))

    addConcentrations(exp, concentrations)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
    expect_true("hasKnownConcentrations" %in% colnames(rowData(exp)))

    expect_true(sum(rowData(exp)$hasKnownConcentrations) == nrow(concentrations))
})



test_that("An experiment can be converted", {
    x <- convertExperiment(exp)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
})


test_that("A message can be sent ", {
    expect_message(message("test"))
})

test_that("an experiment can be converted", {
    expect_true(isValidExperiment(convertExperiment(exp)))
})

test_that("expToCombined", {
    combined <- expToCombined(exp)
    expect_true(is(combined, "data.frame"))
})

test_that("Concentrations can be added", {
    exp <- addConcentrations(exp, read.delim(concentrationsFile))
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(isValidExperiment(exp))
})
