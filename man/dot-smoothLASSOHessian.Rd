% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penaltyFunctions.R
\name{.smoothLASSOHessian}
\alias{.smoothLASSOHessian}
\title{.smoothLASSOHessian}
\usage{
.smoothLASSOHessian(parameters, tuningParameters, penaltyFunctionArguments)
}
\arguments{
\item{parameters}{vector with labeled parameter values}

\item{tuningParameters}{list with field lambda (tuning parameter value)}

\item{penaltyFunctionArguments}{list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)}
}
\value{
Hessian matrix
}
\description{
smoothed version of non-differentiable LASSO Hessian
}
\keyword{internal}
