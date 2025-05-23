% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Concentrations.R
\name{getWeights}
\alias{getWeights}
\title{Get Weights for Weighted Regression}
\usage{
getWeights(yValues)
}
\arguments{
\item{yValues}{A numeric matrix of response values.}
}
\value{
A numeric vector of weights for the regression model.
}
\description{
Computes weights for a weighted regression model based on the
  provided response values. The weights are inversely proportional to the
  values and their squared magnitude.
}
\details{
This function determines the appropriate weights for regression
  by evaluating the response values. If the first column of the response
  matrix contains only zeros, the second column is used instead.
}
