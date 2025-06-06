% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Concentrations.R
\name{modelWithOutliers}
\alias{modelWithOutliers}
\title{Robust Linear Model with Outlier Detection}
\usage{
modelWithOutliers(
  xValues,
  yValues,
  weights,
  forceOrigin = FALSE,
  subtractIntercept = FALSE,
  checkForOutliers = TRUE,
  useWeights = FALSE
)
}
\arguments{
\item{xValues}{A numeric matrix of x-values (e.g., injections or time).}

\item{yValues}{A numeric matrix of y-values (e.g., intensity values).}

\item{weights}{A numeric matrix of weights. Used only if \code{useWeights}
is TRUE.}

\item{forceOrigin}{Logical. If TRUE, fits the model through the origin.}

\item{subtractIntercept}{Logical. If TRUE, subtracts intercept from x-values.}

\item{checkForOutliers}{Logical. If TRUE, detects and removes outliers.}

\item{useWeights}{Logical. If TRUE, computes and applies weights to x-values.}
}
\value{
A list with components:
  \item{slope}{Numeric vector of slopes for each row.}
  \item{intercept}{Numeric vector of intercepts for each row.}
  \item{residuals}{Matrix of residuals for each sample.}
  \item{outliers}{Logical matrix indicating outliers.}
}
\description{
Fits a row-wise linear model between x and y, with options
  for intercept handling, outlier detection using studentized residuals,
  and reweighting.
}
\details{
This function calculates slope and intercept row-wise and
  optionally removes outliers from the first and last samples before
  refitting. Studentized residuals greater than 2 are considered outliers.
}
