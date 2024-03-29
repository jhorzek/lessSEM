% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{.compileTransformations}
\alias{.compileTransformations}
\title{.compileTransformations}
\usage{
.compileTransformations(syntax, parameterLabels, compile = TRUE, notes = NULL)
}
\arguments{
\item{syntax}{string with user defined transformations}

\item{parameterLabels}{names of parameters in the model}

\item{compile}{if set to FALSE, the function will not be compiled -> for visual inspection}

\item{notes}{option to pass a notes to function. All notes of the current
function will be added}
}
\value{
list with parameter names and two Rcpp functions: (1) the transformation function and
(2) a function to create a pointer to the transformation function.
If starting values were defined, these are returned as well.
}
\description{
compile user defined parameter transformations to a
pass to a SEM
}
\keyword{internal}
