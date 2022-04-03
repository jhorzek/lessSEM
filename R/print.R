# https://stackoverflow.com/questions/47871377/print-generic-for-rcpp-class
setClass("Rcpp_SEMCpp")
setMethod("show", "Rcpp_SEMCpp", function (SEM) {
  cat("Internal C++ model representation of aCV4SEM\n")
  cat("Parameters:\n")
  print(aCV4SEM::getParameters(SEM))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", SEM$m2LL))
})