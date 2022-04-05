#' aCV4SEM
#'
#' aCV4SEM implements approximate cross-validation for SEM
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib aCV4SEM
#' @name aCV4SEM

Rcpp::loadModule("SEM_cpp", TRUE)

# https://stackoverflow.com/questions/47871377/print-generic-for-rcpp-class
setMethod("show", "Rcpp_SEMCpp", function (object) {
  cat("Internal C++ model representation of aCV4SEM\n")
  cat("Parameters:\n")
  print(aCV4SEM::getParameters(object))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", object$m2LL))
})

setMethod("logLik", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- -.5*object$m2LL
  class(ll) <- "logLik"
  return(ll)
})
