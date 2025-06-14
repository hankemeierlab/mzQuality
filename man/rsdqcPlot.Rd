% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plots.R
\name{rsdqcPlot}
\alias{rsdqcPlot}
\title{Plot heatmap of RSDQCs with internal standards}
\usage{
rsdqcPlot(exp)
}
\arguments{
\item{exp}{SummarizedExperiment object}
}
\value{
Plotly heatmap of compound - internal standard RSDQCs
}
\description{
This function plots a heatmap of batch-corrected RSDQCs
of all compound - internal standard combinations. This allows for insight
if the picked internal standard is the correct one or if possible another
internal standard should be picked.
}
\examples{
path <- system.file("extdata", "example.tsv", package = "mzQuality")
exp <- buildExperiment(readData(path))

# Run analysis
exp <- doAnalysis(exp)

# Plot RSDQCs
rsdqcPlot(exp)
}
