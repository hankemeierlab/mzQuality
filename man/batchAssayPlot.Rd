% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plots.R
\name{batchAssayPlot}
\alias{batchAssayPlot}
\title{Plot an assay per batch with point distribution}
\usage{
batchAssayPlot(exp, assay = 1, compound = 1)
}
\arguments{
\item{exp}{SummarizedExperiment object}

\item{assay}{Name or index of the assay to plot}

\item{compound}{Compound to plot}
}
\value{
ggplot object of a Batch Assay Plot
}
\description{
Plot an assay per batch with point distribution
}
\examples{
path <- system.file("extdata", "example.tsv", package = "mzQuality")
exp <- buildExperiment(readData(path))

# Create Batch Assay Plot
batchAssayPlot(exp)
}
