#' @title mzQuality Package
#' @description This is an R package for monitoring technical variation in mass
#' spectrometry measurements. It provides various metrics by utilizing pooled
#' quality control samples, blank samples and calibration lines. mzQuality uses
#' Bioconductors' SummarizedExperiment to facilitate implementation with
#' existing workflows. It is fledged with an Rshiny application to support
#' interactive assessment of measurements.
#' @keywords iternal, Software, MassSpectrometry, BatchEffect,
#' QualityControl, Metabolomics
"_PACKAGE"

#' @title Create a timestamped message
#' @description
#' This overrides the base message with an added timestamp
message <- function(message){
    string <- sprintf("%s | %s", lubridate::now(), message)
    base::message(string)
}

pkg.env <- new.env()
pkg.env$rowDataExclude <- c(
    "calModel",
    "concentrationR2",
    "studentizedResiduals",
    "concentrationOutliers",
    "linearRanges",
    "calRatios",
    "compound"
)

pkg.env$colDataExclude <- c(
    "color"
)

pkg.env$assayExclude <- c(
    "ACALRange"
)
