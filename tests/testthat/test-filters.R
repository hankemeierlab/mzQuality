library(testthat)
# Read the example dataset
exp <- readRDS(system.file("data.RDS", package = "mzQuality"))
exp <- doAnalysis(exp, doAll = TRUE, removeOutliers = TRUE)

test_that("filterISTD works", {
    # Should be no internal standards left in the compounds
    exp <- filterISTD(exp, tag = "ISTD")

    # Check if the experiment is still valid
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(isValidExperiment(exp))
})

test_that("filterSST works", {
    # Should be no System Sustainability Tests (SSTs) left in the aliquos
    exp <- filterSST(exp, tag = "SST")

    # Check if the experiment is still valid
    expect_true(is(exp, "SummarizedExperiment"))
    expect_true(isValidExperiment(exp))
})

test_that("filterRSDQC works", {
    rsdqc <- rowData(exp)$rsdqcCorrected
    highConf <- sum(rsdqc <= 15)
    x <- filterRSDQC(exp, min = 0, max = 15)
    expect_equal(nrow(x), highConf)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))

    lowConf <- sum(rsdqc > 15 & rsdqc <= 30)
    x <- filterRSDQC(exp, min = 15, max = 30)
    expect_equal(nrow(x), lowConf)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))

    noConf <- sum(rsdqc > 30)
    x <- filterRSDQC(exp, min = 30, max = Inf)
    expect_equal(nrow(x), noConf)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
})

test_that("filterBackground works", {
    background <- rowData(exp)$backgroundSignal
    highConf <- sum(background <= 0.4)
    x <- filterBackground(exp, min = 0, max = 0.4)
    expect_equal(nrow(x), highConf)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
})

test_that("filterOutliers works", {
    x <- filterOutliers(exp)
    expect_true(is(x, "SummarizedExperiment"))
    expect_true(isValidExperiment(x))
})

