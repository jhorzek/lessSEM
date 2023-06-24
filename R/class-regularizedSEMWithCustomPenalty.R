#' Class for regularized SEM using Rsolnp
#' @slot parameters data.frame with parameter estimates
#' @slot fits data.frame with all fit results
#' @slot parameterLabels character vector with names of all parameters
#' @slot internalOptimization list of elements used internally
#' @slot inputArguments list with elements passed by the user to the general
#' @slot notes internal notes that have come up when fitting the model
#' @keywords internal
setClass(Class = "regularizedSEMWithCustomPenalty",
         representation = representation(
           parameters="data.frame",
           fits="data.frame", 
           parameterLabels = "character",
           internalOptimization = "list",
           inputArguments="list",
           notes = "character"
         )
)

#' summary
#' @param object object of class regularizedSEMWithCustomPenalty
#' @param ... not used
#' @return No return value, just prints estimates
#' @keywords internal
setMethod("summary", "regularizedSEMWithCustomPenalty", function (object, ...) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class regularizedSEMWithCustomPenalty ####\n\n"))
  cat(paste0("- Use coef(", modelName, 
             ") to get the parameter estimates of the model. With coef(", 
             modelName, "lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Information criteria can be computed with AIC(", modelName, 
             ") or BIC(", modelName, 
             ").\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of a regularizedSEMWithCustomPenalty
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @param ... not used
#' @returns data.frame with all parameter estimates
#' @keywords internal
setMethod("coef", "regularizedSEMWithCustomPenalty", function (object, ...) {
  tuningParameters <- object@parameters[, !colnames(object@parameters) %in% object@parameterLabels,drop=FALSE] 
  estimates <- as.matrix(object@parameters[,object@parameterLabels,drop=FALSE])

  if(ncol(object@transformations) != 0){
    transformations <- as.matrix(object@transformations[,
                                                        !colnames(object@transformations) %in% colnames(tuningParameters), 
                                                        drop = FALSE])
  }else{
    transformations <- matrix(nrow = 0, ncol = 0)
  }
  
  coefs <- new("lessSEMCoef")
  coefs@tuningParameters <- tuningParameters
  coefs@estimates <- estimates
  coefs@transformations <- transformations
  
  return(coefs)
})

#' AIC
#' 
#' returns the AIC. Expects penalizedParameterLabels and zeroThreshold
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @param ... Expects penalizedParameterLabels and zeroThreshold. 
#' penalizedParameterLabels: vector with labels of penalized parameters.
#' zeroThreshold: penalized parameters below this threshold will be counted as zeroed.
#' @param k multiplier for number of parameters
#' @returns AIC values
#' @keywords internal
setMethod("AIC", "regularizedSEMWithCustomPenalty", function (object, ..., k) {
  
  dots <- list(...)
  if(!"penalizedParameterLabels" %in% names(dots))
    stop("Expected penalizedParameterLabels")
  if(!"zeroThreshold" %in% names(dots))
    stop("Expected zeroThreshold")
  
  penalizedParameterLabels <- dots$penalizedParameterLabels
  zeroThreshold <- dots$zeroThreshold
  
  fits <- object@fits
  parameters <- object@parameters
  
  fits$AIC <- rep(NA, length(fits$m2LL))
  
  fits$AIC <- fits$m2LL + k*apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  fits$nonZeroParameters <- apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  return(fits)
  
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @param ... Expects penalizedParameterLabels and zeroThreshold. 
#' penalizedParameterLabels: vector with labels of penalized parameters.
#' zeroThreshold: penalized parameters below this threshold will be counted as zeroed.
#' @returns BIC values
#' @keywords internal
setMethod("BIC", "regularizedSEMWithCustomPenalty", function (object, ...) {
  dots <- list(...)
  if(!"penalizedParameterLabels" %in% names(dots))
    stop("Expected penalizedParameterLabels")
  if(!"zeroThreshold" %in% names(dots))
    stop("Expected zeroThreshold")
  
  penalizedParameterLabels <- dots$penalizedParameterLabels
  zeroThreshold <- dots$zeroThreshold
  
  fits <- object@fits
  parameters <- object@parameters
  
  fits$BIC <- rep(NA, length(fits$m2LL))
  
  N <- lavaan::lavInspect(object = object@inputArguments$lavaanModel, what = "nobs")
  
  fits$BIC <- fits$m2LL + log(N)*apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  fits$nonZeroParameters <- apply(parameters, 1, function(x) length(x) - sum(abs(x[penalizedParameterLabels]) < zeroThreshold))
  return(fits)
})
