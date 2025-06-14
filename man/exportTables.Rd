% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Export.R
\name{exportTables}
\alias{exportTables}
\title{Export rowData, colData, and Assays}
\usage{
exportTables(exp, folder)
}
\arguments{
\item{exp}{A SummarizedExperiment object.}

\item{folder}{A folder where the tables will be stored.}
}
\value{
Creates files of rowData, colData, and assays stored in `exp`.
}
\description{
Produces human-readable tables of rowData, colData, and assays
  from a SummarizedExperiment object. The rowData is stored as
  'Compounds.tsv', colData as 'Aliquots.tsv', and assays are stored under
  their respective assay names.
}
\details{
This function creates a folder named "Export" in the specified
  directory and writes the rowData, colData, and assays of the
  SummarizedExperiment object into separate tab-separated files.
}
\examples{
path <- system.file("extdata", "example.tsv", package = "mzQuality")
exp <- buildExperiment(readData(path))

# Perform analysis
exp <- doAnalysis(exp)

# Export tables in the experiment
exportTables(exp[1:10, 1:10], folder = tempdir())
}
