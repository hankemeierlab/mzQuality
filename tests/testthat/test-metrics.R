library(testthat)
suppressPackageStartupMessages(library(mzQuality))
path <- system.file("extdata", "example.tsv", package = "mzQuality")
exp <- buildExperiment(readData(path))
exp <- doAnalysis(exp, removeOutliers = TRUE)

test_that("Outliers are detected correctly", {
    # For each assay, there should be two outliers
    assays <- c("ratio", "area", "ratio_corrected")

    for (a in assays) {
        # Run the outlier detection
        x <- identifyOutliers(exp, assay = a)

        # Check if still is a SummarizedExperiment
        expect_true(is(x, "SummarizedExperiment"))
        expect_true(isValidExperiment(x))
        expect_true("outlier" %in% colnames(colData(x)))
    }
})

test_that("identifyMisInjections", {
    typeName <- "SAMPLE"
    x <- identifyMisInjections(exp, assay = "area_is", type = typeName)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))

    # No misinjections are expected
    expect_true(!all(x[, x$type == typeName]$use))

    # Non-existing sample should return the original experiment
    typeName <- "wrongType"

    x <- identifyMisInjections(exp, assay = "area", type = typeName)
    expect_true(sum(colData(x)$type[x$outlier] == typeName) == 0)
})

test_that("ratioQcSample", {
    exp <- ratioQcSample(exp)
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(isValidExperiment(exp))
    expect_true("ratioQcSample" %in% colnames(rowData(exp)))
})

test_that("doAnalysis is exectued properly", {
    x <- doAnalysis(exp, useWithinBatch = TRUE, removeOutliers = TRUE)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
})

test_that("CarryOver effect can be calculated", {
    x <- carryOverEffect(exp)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
})

test_that("Batch effects are calculated correctly", {
    x <- addBatchCorrection(exp, assay = "ratio")
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
})
