# https://gallery.rcpp.org/articles/custom-printer-exposed-modules/
#' internal representation of SEM in C++
setClass("Rcpp_mgSEM")

#' show
#' 
#' @param object object of class Rcpp_mgSEM
setMethod("show", "Rcpp_mgSEM", function (object) {
  cat("Internal C++ model representation of lessSEM\n")
  cat("Parameters:\n")
  print(.getParameters(object))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", object$m2LL))
})

#' logLik
#' 
#' @param object object of class Rcpp_mgSEM
#' @importMethodsFrom stats logLik
#' @returns log-likelihood of the model
setMethod("logLik", "Rcpp_mgSEM", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  N <- object$sampleSize
  numberOfParameters <- length(.getParameters(object))
  
  ll <- new("logLikelihood",
            logLik = -.5*object$m2LL,
            nParameters = numberOfParameters,
            N = N)
  return(ll)
})

#' coef
#' 
#' @param object object of class Rcpp_mgSEM
#' @param ... not used
#' @importMethodsFrom stats coef
#' @returns all coefficients of the model in transformed form
setMethod("coef", "Rcpp_mgSEM", function (object, ...) {
  return(.getParameters(object, raw = FALSE))
})

#' AIC
#' 
#' @param object object of class Rcpp_mgSEM
#' @importMethodsFrom stats AIC
#' @returns AIC values
setMethod("AIC", "Rcpp_mgSEM", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  AICis <- -2*ll@logLik + 2*ll@nParameters
  return(AICis)
})

#' BIC
#' 
#' @param object object of class Rcpp_mgSEM
#' @importMethodsFrom stats BIC
#' @returns BIC values
setMethod("BIC", "Rcpp_mgSEM", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  BICis <- -2*ll@logLik + log(ll@N)*ll@nParameters
  return(BICis)
})
