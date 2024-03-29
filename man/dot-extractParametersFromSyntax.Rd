% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{.extractParametersFromSyntax}
\alias{.extractParametersFromSyntax}
\title{.extractParametersFromSyntax}
\usage{
.extractParametersFromSyntax(syntax, parameterLabels)
}
\arguments{
\item{syntax}{syntax for parameter transformations}

\item{parameterLabels}{names of parameters in the model}
}
\value{
vector with names of parameters used in the syntax and vector with
boolean indicating if parameter is transformation result
}
\description{
extract the names of the parameters in a syntax
}
\keyword{internal}
