% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plots.R
\name{rsdPlot}
\alias{rsdPlot}
\title{Plot the Relative Standard Deviation}
\usage{
rsdPlot(
  exp,
  compound = 1,
  assay = "ratio_corrected",
  qc = "SQC",
  number = nrow(exp)
)
}
\arguments{
\item{exp}{SummarizedExperiment object}

\item{compound}{Which compound to plot. Defaults to the first compound.}

\item{assay}{Which assay to use for calculating the RSDQC}

\item{qc}{type of QC to plot}

\item{number}{Optional Show the n number of worst RSD}
}
\value{
ggplot object of a scatter- / lineplot of RSDs per compound
and batch
}
\description{
This plot shows the RSD of compounds across batches. This can
aid in identifying which compounds have a high variation in the given assay
in certain batches.
}
\examples{
path <- system.file("extdata", "example.tsv", package = "mzQuality")
exp <- buildExperiment(readData(path))

# Do analysis
exp <- doAnalysis(exp)

# Create scatterplot with Relative Standard Deviation values
rsdPlot(exp, assay = "ratio", number = 10)
}
