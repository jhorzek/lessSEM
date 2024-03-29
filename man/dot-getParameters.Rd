% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{.getParameters}
\alias{.getParameters}
\title{.getParameters}
\usage{
.getParameters(SEM, raw = FALSE, transformations = FALSE)
}
\arguments{
\item{SEM}{model of class Rcpp_SEMCpp. Models of this class}

\item{raw}{controls if the parameter are returned in raw format or transformed}

\item{transformations}{should transformed parameters be included?}
}
\value{
labeled vector with parameter values
}
\description{
returns the parameters of the internal model representation.
}
\keyword{internal}
