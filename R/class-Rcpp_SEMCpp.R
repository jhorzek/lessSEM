# https://gallery.rcpp.org/articles/custom-printer-exposed-modules/
#' internal representation of SEM in C++
#' @keywords internal
setClass("Rcpp_SEMCpp")

#' show
#' 
#' @param object object of class Rcpp_SEMCpp
#' @return No return value, just prints estimates
setMethod("show", "Rcpp_SEMCpp", function (object) {
  cat("Internal C++ model representation of lessSEM\n")
  cat("Parameters:\n")
  print(.getParameters(object))
  cat("\n")
  cat(paste0("Objective value: ", object$objectiveValue))
})

#' logLik
#' 
#' @param object object of class Rcpp_SEMCpp
#' @param ... not used
#' @returns log-likelihood of the model
setMethod("logLik", "Rcpp_SEMCpp", function (object, ...) {
  
  if(object$getEstimator() != "fiml")
    stop("Likelihood can only be computed for fiml objective.")
  
  if(!object$wasFit){
    object$fit()
  }
  N <- nrow(object$rawData)
  numberOfParameters <- length(.getParameters(object))
  
  ll <- new("logLikelihood",
            logLik = -.5*object$objectiveValue,
            nParameters = numberOfParameters,
            N = N)
  return(ll)
})

#' coef
#' 
#' @param object object of class Rcpp_SEMCpp
#' @param ... not used
#' @returns all coefficients of the model in transformed form
setMethod("coef", "Rcpp_SEMCpp", function (object, ...) {
  return(.getParameters(object, raw = FALSE))
})

#' AIC
#' 
#' @param object object of class Rcpp_SEMCpp
#' @param ... not used
#' @param k multiplier for number of parameters
#' @returns AIC values
setMethod("AIC", "Rcpp_SEMCpp", function (object, ..., k = 2) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  AICis <- -2*ll@logLik + k*ll@nParameters
  return(AICis)
})

#' BIC
#' 
#' @param object object of class Rcpp_SEMCpp
#' @param ... not used
#' @returns BIC values
setMethod("BIC", "Rcpp_SEMCpp", function (object, ...) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  BICis <- -2*ll@logLik + log(ll@N)*ll@nParameters
  return(BICis)
})
