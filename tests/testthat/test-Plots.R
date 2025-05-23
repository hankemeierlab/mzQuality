library(testthat)
suppressPackageStartupMessages(library(mzQuality))
# Read the example dataset
exp <- readRDS(system.file("extdata/data.RDS", package = "mzQuality"))
exp <- doAnalysis(exp, useWithinBatch = F)

test_that("pca plot can be made", {
    plot <- pcaPlot(exp)
    expect_true(is(plot, "ggplot"))
})

test_that("Violin Plot can be made", {
    plot <- suppressWarnings(violinPlot(exp))
    expect_true(is(plot, "ggplot"))
})

test_that("Compound plot can be made", {
    plot <- compoundPlot(exp)
    expect_true(is(plot, "ggplot"))

    plot <- compoundPlot(exp, assay = "area")
    expect_true(is(plot, "ggplot"))

    plot <- compoundPlot(exp, assay = "area", addInternalStandards = TRUE)
    expect_true(is(plot, "ggplot"))

})

test_that("Plots can be facetted", {
    plot <- compoundPlot(exp, addText = TRUE)

    p <- facetPlot(plot, shareX = TRUE, shareY = TRUE)
    expect_true(is(p, "ggplot"))

    p <- facetPlot(plot, shareX = FALSE, shareY = FALSE)
    expect_true(is(p, "ggplot"))

    p <- facetPlot(plot, shareX = FALSE, shareY = TRUE)
    expect_true(is(p, "ggplot"))

    p <- facetPlot(plot, shareX = TRUE, shareY = FALSE)
    expect_true(is(p, "ggplot"))

})

test_that("Aliquot Plot can be made", {
    plot <- aliquotPlot(exp, "ratio")
    expect_true(is(plot, "ggplot"))
})

test_that("Concentration Plot can be made", {
    plot <- concentrationPlot(exp)
    expect_true(is(plot, "ggplot"))
})

test_that("Labels can be added to a PCA plot", {
    plot <- pcaPlot(exp, addLabels = TRUE)
    expect_true(is(plot, "ggplot"))
})

test_that("rsdqcPlot can be made", {
    plot <- rsdqcPlot(exp)
    expect_true(is(plot, "plotly"))
})

test_that("heatmapPlot can be made", {
    plot <- heatmapPlot(exp)
    expect_true(is(plot, "plotly"))
})


test_that("rsdPlot can be made", {
    plot <- rsdPlot(exp)
    expect_true(is(plot, "ggplot"))
})

test_that("batchAssayPlot can be made", {
    plot <- batchAssayPlot(exp)
    expect_true(is(plot, "ggplot"))
})


