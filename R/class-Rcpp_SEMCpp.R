# https://gallery.rcpp.org/articles/custom-printer-exposed-modules/
#' internal representation of SEM in C++
setClass("Rcpp_SEMCpp")

#' show
#' 
#' @param object object of class Rcpp_SEMCpp
setMethod("show", "Rcpp_SEMCpp", function (object) {
  cat("Internal C++ model representation of lessSEM\n")
  cat("Parameters:\n")
  print(.getParameters(object))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", object$m2LL))
})

#' logLik
#' 
#' @param object object of class Rcpp_SEMCpp
#' @importMethodsFrom stats logLik
#' @returns log-likelihood of the model
setMethod("logLik", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  N <- nrow(object$rawData)
  numberOfParameters <- length(.getParameters(object))
  
  ll <- new("logLikelihood",
            logLik = -.5*object$m2LL,
            nParameters = numberOfParameters,
            N = N)
  return(ll)
})

#' coef
#' 
#' @param object object of class Rcpp_SEMCpp
#' @param ... not used
#' @importMethodsFrom stats coef
#' @returns all coefficients of the model in transformed form
setMethod("coef", "Rcpp_SEMCpp", function (object, ...) {
  return(.getParameters(object, raw = FALSE))
})

#' AIC
#' 
#' @param object object of class Rcpp_SEMCpp
#' @importMethodsFrom stats AIC
#' @returns AIC values
setMethod("AIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  AICis <- -2*ll@logLik + 2*ll@nParameters
  return(AICis)
})

#' BIC
#' 
#' @param object object of class Rcpp_SEMCpp
#' @importMethodsFrom stats BIC
#' @returns BIC values
setMethod("BIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  BICis <- -2*ll@logLik + log(ll@N)*ll@nParameters
  return(BICis)
})
