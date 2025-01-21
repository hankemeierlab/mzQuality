library(testthat)
suppressPackageStartupMessages(library(mzQuality2))
# Read the example dataset
exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))

test_that("Outliers are detected correctly", {
    # For each assay, there should be two outliers
    assays <- c("ratio", "area", "ratio_corrected")

    for (a in assays) {
        # Run the outlier detection
        x <- identifyOutliers(exp, assay = a)

        # Check if still is a SummarizedExperiment
        expect_true(is(x, "SummarizedExperiment"))
        expect_true(validateExperiment(x))
        expect_true("outlier" %in% colnames(colData(x)))

        # In this dataset there should be 2 outliers:
        expect_true(is.logical(colData(x)$outlier))
        expect_true(sum(colData(x)$outlier) == 2)
    }
})

test_that("identifyMisInjections", {
    typeName <- "SAMPLE"
    x <- identifyMisInjections(exp, assay = "area", type = typeName)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(validateExperiment(x))

    # No misinjections are expected
    expect_true(all(x[, x$type == typeName]$use))

    # Non-existing sample should return the original experiment
    typeName = "wrongType"
    x <- identifyMisInjections(exp, assay = "area", type = typeName)
    expect_equal(x, exp)
})

test_that("ratioQcSample", {

    exp <- ratioQcSample(exp)
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(validateExperiment(exp))
    expect_true("ratioQcSample" %in% colnames(rowData(exp)))
})



# test_that("typePresence", {
#     x <- typePresence(exp, area = "area", type = "SQC")
#     expect_true(all(x == 1))
#     x <- typePresence(exp, area = "area", type = "SAMPLE")
#     expect_true(all(x == 1))
# })

test_that("doAnalysis is exectued properly", {
    x <- doAnalysis(exp, doAll = TRUE, useWithinBatch = TRUE, removeOutliers = TRUE)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(validateExperiment(x))
})

test_that("CarryOver effect can be calculated", {
    x <- carryOverEffect(exp, method = "previous")
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(validateExperiment(x))
})
