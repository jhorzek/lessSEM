% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initializeSEM.R
\name{.multiGroupSEMFromLavaan}
\alias{.multiGroupSEMFromLavaan}
\title{.multiGroupSEMFromLavaan}
\usage{
.multiGroupSEMFromLavaan(
  lavaanModels,
  whichPars = "est",
  fit = TRUE,
  addMeans = TRUE,
  transformations = NULL,
  transformationList = list(),
  transformationGradientStepSize = 1e-06
)
}
\arguments{
\item{lavaanModels}{vector with lavaan models}

\item{whichPars}{which parameters should be used to initialize the model. If set to "est", the parameters will be set to the
estimated parameters of the lavaan model. If set to "start", the starting values of lavaan will be used. The latter can be useful if parameters are to
be optimized afterwards as setting the parameters to "est" may result in the model getting stuck in a local minimum.}

\item{fit}{should the model be fitted}

\item{addMeans}{If lavaanModel has meanstructure = FALSE, addMeans = TRUE will add a mean structure. FALSE will set the means of the observed variables to the average}

\item{transformations}{string with transformations}

\item{transformationList}{list for transformations}

\item{transformationGradientStepSize}{step size used to compute the gradients of the
transformations}
}
\value{
Object of class Rcpp_mgSEMCpp
}
\description{
internal function. Translates a vector of objects of class lavaan to the
internal model representation.
}
\keyword{internal}
