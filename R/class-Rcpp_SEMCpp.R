# https://gallery.rcpp.org/articles/custom-printer-exposed-modules/
setClass("Rcpp_SEMCpp")
setMethod("show", "Rcpp_SEMCpp", function (object) {
  cat("Internal C++ model representation of aCV4SEM\n")
  cat("Parameters:\n")
  print(aCV4SEM:::getParameters(object))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", object$m2LL))
})

setMethod("logLik", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  N <- nrow(object$rawData)
  numberOfParameters <- length(aCV4SEM:::getParameters(object))
  
  ll <- new("logLikelihood",
            logLik = -.5*object$m2LL,
            nParameters = numberOfParameters,
            N = N)
  return(ll)
})

setMethod("coef", "Rcpp_SEMCpp", function (object) {
  return(aCV4SEM:::getParameters(object, raw = FALSE))
})

setMethod("AIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  AICis <- -2*ll@logLik + 2*ll@nParameters
  return(AICis)
})

setMethod("BIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  BICis <- -2*ll@logLik + log(ll@N)*ll@nParameters
  return(BICis)
})
