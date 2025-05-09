% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Metrics.R
\name{typePresence}
\alias{typePresence}
\title{Calculate the presence of a specific type in the experiment}
\usage{
typePresence(exp, area = metadata(exp)$primary, type = metadata(exp)$QC)
}
\arguments{
\item{exp}{A SummarizedExperiment object containing the experimental data.}

\item{area}{A string specifying the assay to use for the calculation.
Defaults to the `primary` metadata of the experiment.}

\item{type}{A string specifying the type of samples to analyze. Defaults
to the `QC` metadata of the experiment.}
}
\value{
A SummarizedExperiment object with the calculated presence added
    to the `rowData`.
}
\description{
Computes the proportion of non-missing values for a specific
    type (e.g., QC) in the assay data of a SummarizedExperiment object.
}
\details{
This function calculates the proportion of non-missing values for
    a specified type (e.g., QC) in the assay data. The result is stored in
    the `rowData` of the experiment under a column named
    `<type>Presence`.
}
