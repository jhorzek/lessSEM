% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penaltyFunctions.R
\name{.ridgeGradient}
\alias{.ridgeGradient}
\title{.ridgeGradient}
\usage{
.ridgeGradient(parameters, tuningParameters, penaltyFunctionArguments)
}
\arguments{
\item{parameters}{vector with labeled parameter values}

\item{tuningParameters}{list with field lambda (tuning parameter value)}

\item{penaltyFunctionArguments}{list with field regularizedParameterLabels (labels of regularized parameters)}
}
\value{
gradient values
}
\description{
ridge gradient function
}
\keyword{internal}
