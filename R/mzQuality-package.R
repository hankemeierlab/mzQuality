#' @title mzQuality Package
#' @description This is an R package for monitoring variation in mass
#' spectrometry measurements. It provides various metrics by utilizing pooled
#' quality control samples, blank samples and calibration lines. mzQuality uses
#' Bioconductors' SummarizedExperiment to facilitate implementation with
#' existing workflows.
#' @keywords iternal, Software, MassSpectrometry, BatchEffect,
#' QualityControl, Metabolomics
"_PACKAGE"

pkg.env <- new.env()
pkg.env$rowDataExclude <- c(
    "concentrationR2",
    "studentizedResiduals",
    "linearRanges",
    "calRatios",
    "compound"
)

pkg.env$assayExclude <- c(
    "CALRange",
    "ACALRange"
)
