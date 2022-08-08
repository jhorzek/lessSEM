# https://gallery.rcpp.org/articles/custom-printer-exposed-modules/
#' internal representation of SEM in C++
#' 
setClass("Rcpp_SEMCpp")

#' show
#' 
#' @param object object of class Rcpp_SEMCpp
#' @export
setMethod("show", "Rcpp_SEMCpp", function (object) {
  cat("Internal C++ model representation of lessSEM\n")
  cat("Parameters:\n")
  print(lessSEM:::.getParameters(object))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", object$m2LL))
})

#' logLik
#' 
#' @param object object of class Rcpp_SEMCpp
#' @returns log-likelihood of the model
#' @export
setMethod("logLik", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  N <- nrow(object$rawData)
  numberOfParameters <- length(lessSEM:::.getParameters(object))
  
  ll <- new("logLikelihood",
            logLik = -.5*object$m2LL,
            nParameters = numberOfParameters,
            N = N)
  return(ll)
})

#' coef
#' 
#' @param object object of class Rcpp_SEMCpp
#' @returns all coefficients of the model in transformed form
#' @export
setMethod("coef", "Rcpp_SEMCpp", function (object) {
  return(lessSEM:::.getParameters(object, raw = FALSE))
})

#' AIC
#' 
#' @param object object of class Rcpp_SEMCpp
#' @returns AIC values
#' @export
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
#' @returns BIC values
#' @export
setMethod("BIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  BICis <- -2*ll@logLik + log(ll@N)*ll@nParameters
  return(BICis)
})
