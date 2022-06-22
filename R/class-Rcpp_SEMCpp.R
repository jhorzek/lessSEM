# https://gallery.rcpp.org/articles/custom-printer-exposed-modules/
setClass("Rcpp_SEMCpp")

#' @export
setMethod("show", "Rcpp_SEMCpp", function (object) {
  cat("Internal C++ model representation of linr\n")
  cat("Parameters:\n")
  print(linr:::getParameters(object))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", object$m2LL))
})

#' @export
setMethod("logLik", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  N <- nrow(object$rawData)
  numberOfParameters <- length(linr:::getParameters(object))
  
  ll <- new("logLikelihood",
            logLik = -.5*object$m2LL,
            nParameters = numberOfParameters,
            N = N)
  return(ll)
})

#' @export
setMethod("coef", "Rcpp_SEMCpp", function (object) {
  return(linr:::getParameters(object, raw = FALSE))
})

#' @export
setMethod("AIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  AICis <- -2*ll@logLik + 2*ll@nParameters
  return(AICis)
})

#' @export
setMethod("BIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  BICis <- -2*ll@logLik + log(ll@N)*ll@nParameters
  return(BICis)
})
