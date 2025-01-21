library(testthat)
suppressPackageStartupMessages(library(mzQuality2))
# Read the example dataset
exp <- readRDS(system.file("data.RDS", package = "mzQuality2"))

test_that("pca plot can be made", {
    plot <- pcaPlotNew(exp)
    expect_true(is(plot , "ggplot"))
})

test_that("Violin Plot can be made", {
    plot <- violinPlotNew(exp)
    expect_true(is(plot , "ggplot"))
})

test_that("Compound plot can be made", {
    plot <- compoundPlotNew(exp)
    expect_true(is(plot , "ggplot"))

    plot <- compoundPlotNew(exp, assay = "area")
    expect_true(is(plot , "ggplot"))
})

test_that("Plots can be facetted", {
    plot <- compoundPlotNew(exp) %>%
        facetPlot()
    expect_true(is(plot , "ggplot"))
})

test_that("Aliquot Plot can be made", {
    plot <- aliquotPlotNew(exp, "ratio")
    expect_true(is(plot , "ggplot"))

})

test_that("Concentration Plot can be made", {
    plot <- concentrationPlotNew(exp)
    expect_true(is(plot , "ggplot"))
})

test_that('Labels can be added to a PCA plot', {
    plot <- pcaPlotNew(exp) %>%
        addLabels()
    expect_true(is(plot , "ggplot"))

})

test_that('rsdqcPlot can be made', {
    plot <- rsdqcPlot(exp)
    expect_true(is(plot , "plotly"))
})

test_that('heatmapPlot can be made', {
    plot <- heatmapPlot(exp)
    expect_true(is(plot , "plotly"))
})

test_that('calibrationPlot can be made', {
    plot <- calibrationPlot(exp, assay = "concentration")
    expect_true(is(plot , "ggplot"))
})


test_that("rsdPlot can be made", {
    plot <- rsdPlot(exp)
    expect_true(is(plot , "ggplot"))
})

test_that("batchAssayPlot can be made", {
    plot <- batchAssayPlot(exp)
    expect_true(is(plot , "ggplot"))
})

test_that("batchCorrectionPlot can be made", {
    plot <- batchCorrectionPlot(exp)
    expect_true(is(plot , "ggplot"))
})
