% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Convert.R
\name{expToCombined}
\alias{expToCombined}
\title{Convert a SummarizedExperiment to a data frame in the long format}
\usage{
expToCombined(exp, compoundColumn = "compound", aliquotColumn = "aliquot")
}
\arguments{
\item{exp}{a SummarizedExperiment object}

\item{compoundColumn}{Columnname in the resulting dataframe for compounds,
defaults to 'Compound'.}

\item{aliquotColumn}{Columnname in the resulting dataframe for aliquots,
defaults to 'Aliquot'.}
}
\value{
a data.frame in a generic, human-readable format
}
\description{
A SummarizedExperiment is great for analyzing data, but not
useful for storage in human-readable format. To facilitate in this need,
this function can be used to convert a SummarizedExperiment to a dataframe
in the long format. Note that metadata won't be exported.
}
\examples{
path <- system.file("extdata/example.tsv", package = "mzQuality")
data <- readData(path)

# Construct experiment
exp <- buildExperiment(
    data,
    compoundColumn = "compound",
    aliquotColumn = "aliquot",
    primaryAssay = "area",
    secondaryAssay = "area_is"
)

# Convert to a table format
expToCombined(exp, compoundColumn = "Compound", aliquotColumn = "Aliquot")
}
