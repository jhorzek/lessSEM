% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{istaLSPSEM}
\alias{istaLSPSEM}
\title{lsp optimization with ista}
\value{
a list with fit results
}
\description{
Object for lsp optimization with
ista optimizer
}
\section{Fields}{

\describe{
\item{\code{new}}{creates a new object. Requires (1) a vector with weights for each
parameter and (2) a list with control elements}

\item{\code{optimize}}{optimize the model. Expects a vector with starting values,
a SEM of type SEM_Cpp, a theta and a lambda value.}
}}

