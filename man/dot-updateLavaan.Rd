% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lavaanHelperFunctions.R
\name{.updateLavaan}
\alias{.updateLavaan}
\title{.updateLavaan}
\usage{
.updateLavaan(lavaanModel, key, value)
}
\arguments{
\item{lavaanModel}{fitted lavaan model}

\item{key}{label of the element that should be updated}

\item{value}{new value for the updated element}
}
\value{
lavaan model
}
\description{
updates a lavaan model. lavaan has an update function that does exactly that,
but it seems to not work with testthat. This is an attempt to hack around the
issue...
}
