#' Class for regularized SEM using Rsolnp
#' @slot parameters data.frame with parameter estimates
#' @slot fits data.frame with all fit results
#' @slot parameterLabels character vector with names of all parameters
#' @slot internalOptimization list of elements used internally
#' @slot inputArguments list with elements passed by the user to the general
setClass(Class = "regularizedSEMWithCustomPenalty",
         representation = representation(
           parameters="data.frame",
           fits="data.frame", 
           parameterLabels = "character",
           internalOptimization = "list",
           inputArguments="list"
         )
)

#' summary
#' @param object object of class regularizedSEMWithCustomPenalty
#' @export
setMethod("summary", "regularizedSEMWithCustomPenalty", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class regularizedSEMWithCustomPenalty ####\n\n"))
  cat(paste0("- Use coef(", modelName, 
             ") to get the parameter estimates of the model. With coef(", 
             modelName, "lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Information criteria can be compute with AIC(", modelName, 
             ") or BIC(", modelName, 
             ").\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of a regularizedSEMWithCustomPenalty
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @returns data.frame with all parameter estimates
#' @export
setMethod("coef", "regularizedSEMWithCustomPenalty", function (object, ...) {
  pars <- object@parameters
  return(pars)
})

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @param penalizedParameterLabels vector with labels of penalized parameters
#' @param zeroThreshold penalized parameters below this threshold will be counted as zeroed
#' @returns AIC values
#' @export
setMethod("AIC", "regularizedSEMWithCustomPenalty", function (object, penalizedParameterLabels, zeroThreshold) {
  
  fits <- object@fits
  parameters <- object@parameters
  
  fits$AIC <- rep(NA, length(fits$m2LL))
  
  fits$AIC <- fits$m2LL + 2*apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  fits$nonZeroParameters <- apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  return(fits)
  
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @param penalizedParameterLabels vector with labels of penalized parameters
#' @param zeroThreshold penalized parameters below this threshold will be counted as zeroed
#' @returns BIC values
#' @export
setMethod("BIC", "regularizedSEMWithCustomPenalty", function (object, penalizedParameterLabels, zeroThreshold) {
  fits <- object@fits
  parameters <- object@parameters
  
  fits$BIC <- rep(NA, length(fits$m2LL))
  
  N <- lavaan::lavInspect(object = object@inputArguments$lavaanModel, what = "nobs")
  
  fits$BIC <- fits$m2LL + log(N)*apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  fits$nonZeroParameters <- apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  return(fits)
})
