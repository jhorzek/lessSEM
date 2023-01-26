#' Class for regularized SEM
#' @slot penalty penalty used (e.g., "lasso")
#' @slot parameters data.frame with parameter estimates
#' @slot fits data.frame with all fit results
#' @slot parameterLabels character vector with names of all parameters
#' @slot weights vector with weights given to each of the parameters in the penalty
#' @slot regularized character vector with names of regularized parameters
#' @slot transformations if the model has transformations, the transformed parameters 
#' are returned
#' @slot internalOptimization list of elements used internally
#' @slot inputArguments list with elements passed by the user to the general
setClass(Class = "regularizedSEMMixedPenalty",
         representation = representation(
           penalty = "character",
           tuningParameterConfigurations = "list",
           parameters="data.frame",
           fits="data.frame", 
           parameterLabels = "character",
           weights = "numeric",
           regularized = "character",
           transformations = "data.frame",
           internalOptimization="list", 
           inputArguments="list"
         )
)

#' show
#' @param object object of class regularizedSEM
#' @export
setMethod("show", "regularizedSEMMixedPenalty", function (object) {
  #modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class regularizedSEMMixedPenalty with mixed penalty penalty ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat(paste0("- Use coef(object) to get the parameter estimates of the model.\n\n"))
  cat(paste0("- Information criteria can be computed with AIC(object) or BIC(object).\n\n"))
  cat("################################################\n")
})

#' summary
#' @param object object of class regularizedSEMMixedPenalty
#' @export
setMethod("summary", "regularizedSEMMixedPenalty", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class regularizedSEMMixedPenalty with mixed penalty penalty ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat(paste0("- Use coef(", modelName, 
             ") to get the parameter estimates of the model.\n\n"))
  cat(paste0("- Information criteria can be computed with AIC(", modelName, 
             ") or BIC(", modelName, 
             ").\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of a regularizedSEMMixedPenalty
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @param ... criterion can be one of: "AIC", "BIC". If set to NULL, all parameters will be returned
#' @returns parameters of the model as data.frame
#' @export
setMethod("coef", "regularizedSEMMixedPenalty", function (object, ...) {
  dotdotdot <- list(...)
  if("criterion" %in% names(dotdotdot)){
    criterion <- dotdotdot$criterion
  }else{
    criterion <- NULL
  }
  
  if(!is.null(criterion) && criterion %in% c("AIC", "BIC")){
    if(length(unique(object@fits$nonZeroParameters)) == 1) 
      stop("Selection by criterion currently only supported for sparsity inducing penalties. Either none of your parameters was zeroed or the penalty used does not induce sparsity.")
    if(criterion == "AIC"){
      AICs <- AIC(object)
      bestAIC <- which(AICs$AIC == min(AICs$AIC))[1]
      return(object@parameters[bestAIC,]) 
    }
    
    if(criterion == "BIC"){
      BICs <- BIC(object)
      bestBIC <- which(BICs$BIC == min(BICs$BIC))[1]
      return(object@parameters[bestBIC,])
    }
    
  }
  
  pars <- object@parameters
  return(pars)
})

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @returns AIC values
#' @export
setMethod("AIC", "regularizedSEMMixedPenalty", function (object) {
  
  fits <- object@fits
  fits$AIC <- fits$m2LL + 2*fits$nonZeroParameters
  
  return(fits)
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @returns BIC values
#' @export
setMethod("BIC", "regularizedSEMMixedPenalty", function (object) {
  N <- object@internalOptimization$N
  fits <- object@fits
  
  fits <- object@fits
  fits$BIC <- fits$m2LL + log(N)*fits$nonZeroParameters
  
  return(fits)
  
})