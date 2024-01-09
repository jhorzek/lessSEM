#' Class for regularized SEM
#' @slot penalty penalty used (e.g., "lasso")
#' @slot tuningParameterConfigurations list with settings for the lambda, theta, 
#' and alpha tuning parameters.
#' @slot parameters data.frame with parameter estimates
#' @slot fits data.frame with all fit results
#' @slot parameterLabels character vector with names of all parameters
#' @slot weights vector with weights given to each of the parameters in the penalty
#' @slot regularized character vector with names of regularized parameters
#' @slot transformations if the model has transformations, the transformed parameters 
#' are returned
#' @slot internalOptimization list of elements used internally
#' @slot inputArguments list with elements passed by the user to the general
#' @slot notes internal notes that have come up when fitting the model
#' @export
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
           inputArguments="list",
           notes = "character"
         )
)

#' show
#' @param object object of class regularizedSEM
#' @return No return value, just prints estimates
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
#' @param ... not used
#' @return No return value, just prints estimates
#' @export
setMethod("summary", "regularizedSEMMixedPenalty", function (object, ...) {
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
  
  
  tuningParameters <- object@parameters[, !colnames(object@parameters) %in% object@parameterLabels,drop=FALSE] 
  estimates <- as.matrix(object@parameters[,object@parameterLabels,drop=FALSE])
  
  if(ncol(object@transformations) != 0){
    transformations <- as.matrix(object@transformations[,
                                                        !colnames(object@transformations) %in% colnames(tuningParameters), 
                                                        drop = FALSE])
  }else{
    transformations <- matrix(nrow = 0, ncol = 0)
  }
  
  if(!is.null(criterion)){
    
    fits <- fitIndices(object)
    
    if(!criterion %in% colnames(fits))
      stop("Could not find", criterion, "in fitIndices(object).")
    
    if(length(unique(object@fits$nonZeroParameters)) == 1) 
      stop("Selection by criterion currently only supported for sparsity inducing penalties. Either none of your parameters was zeroed or the penalty used does not induce sparsity.")
    
    bestFit <- which(fits[,criterion] == min(fits[,criterion]))[1]
    
    coefs <- new("lessSEMCoef")
    coefs@tuningParameters <- tuningParameters[bestFit,,drop = FALSE]
    coefs@estimates <- estimates[bestFit,,drop = FALSE]
    if(ncol(object@transformations) != 0){
      coefs@transformations <- transformations[bestFit,,drop = FALSE]
    }else{
      coefs@transformations <- transformations
    }
    return(coefs)
  }
  
  coefs <- new("lessSEMCoef")
  coefs@tuningParameters <- tuningParameters
  coefs@estimates <- estimates
  coefs@transformations <- transformations
  
  return(coefs)
})

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @param ... not used
#' @param k multiplier for number of parameters
#' @returns AIC values
#' @export
setMethod("AIC", "regularizedSEMMixedPenalty", function (object, ..., k = 2) {
  
  fits <- object@fits
  fits$AIC <- fits$m2LL + k*fits$nonZeroParameters
  
  return(fits)
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @param ... not used
#' @returns BIC values
#' @export
setMethod("BIC", "regularizedSEMMixedPenalty", function (object, ...) {
  N <- object@internalOptimization$N
  fits <- object@fits
  
  fits <- object@fits
  fits$BIC <- fits$m2LL + log(N)*fits$nonZeroParameters
  
  return(fits)
  
})

#' getTuningParameterConfiguration
#' 
#' Returns the lambda, theta, and alpha values for the tuning parameters
#' of a regularized SEM with mixed penalty.
#' @param regularizedSEMMixedPenalty object of type regularizedSEMMixedPenalty (see ?mixedPenalty)
#' @param tuningParameterConfiguration integer indicating which tuningParameterConfiguration 
#' should be extracted (e.g., 1). See the entry in the row tuningParameterConfiguration of
#' regularizedSEMMixedPenalty@fits and regularizedSEMMixedPenalty@parameters.
#' @return data frame with penalty and tuning parameter settings
#' @examples 
#' library(lessSEM)
#' 
#' # Identical to regsem, lessSEM builds on the lavaan
#' # package for model specification. The first step
#' # therefore is to implement the model in lavaan.
#' 
#' dataset <- simulateExampleData()
#' 
#' lavaanSyntax <- "
#' f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
#'      l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
#'      l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
#' f ~~ 1*f
#' "
#' 
#' lavaanModel <- lavaan::sem(lavaanSyntax,
#'                            data = dataset,
#'                            meanstructure = TRUE,
#'                            std.lv = TRUE)
#' 
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addLsp(regularized = paste0("l", 11:15), 
#'          lambdas = seq(0,1,.1),
#'          thetas = 2.3) |>
#'   # fit the model:
#'   fit()
#' 
#' getTuningParameterConfiguration(regularizedSEMMixedPenalty = regularized, 
#'                                 tuningParameterConfiguration = 2)
#' @export
getTuningParameterConfiguration <- function(regularizedSEMMixedPenalty, 
                                            tuningParameterConfiguration){
  if(!is(regularizedSEMMixedPenalty, "regularizedSEMMixedPenalty"))
    stop("regularizedSEMMixedPenalty must be of class 'regularizedSEMMixedPenalty'")
  if(length(tuningParameterConfiguration) != 1)
    stop("tuningParameterConfiguration must be a single value (e.g., 1)")
  if(!tuningParameterConfiguration %in% 1:nrow(regularizedSEMMixedPenalty@tuningParameterConfigurations$lambda)){
    stop("tuningParameterConfiguration must be an integer in [1, ", nrow(regularizedSEMMixedPenalty@tuningParameterConfigurations$lambda), "]")
  }
  
  if(regularizedSEMMixedPenalty@inputArguments$method == "glmnet")
    return(
      data.frame(
        parameter = regularizedSEMMixedPenalty@parameterLabels,
        penalty = regularizedSEMMixedPenalty@penalty[regularizedSEMMixedPenalty@parameterLabels],
        lambda = regularizedSEMMixedPenalty@tuningParameterConfigurations$lambda[tuningParameterConfiguration,regularizedSEMMixedPenalty@parameterLabels],
        alpha = regularizedSEMMixedPenalty@tuningParameterConfigurations$alpha[tuningParameterConfiguration,regularizedSEMMixedPenalty@parameterLabels]
      )
    )
  
  return(
    data.frame(
      parameter = regularizedSEMMixedPenalty@parameterLabels,
      penalty = regularizedSEMMixedPenalty@penalty[regularizedSEMMixedPenalty@parameterLabels],
      lambda = regularizedSEMMixedPenalty@tuningParameterConfigurations$lambda[tuningParameterConfiguration,regularizedSEMMixedPenalty@parameterLabels],
      theta = regularizedSEMMixedPenalty@tuningParameterConfigurations$theta[tuningParameterConfiguration,regularizedSEMMixedPenalty@parameterLabels],
      alpha = regularizedSEMMixedPenalty@tuningParameterConfigurations$alpha[tuningParameterConfiguration,regularizedSEMMixedPenalty@parameterLabels]
    )
  )
}

#' fitIndices
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @return returns a data.frame with fit indices
#' @export
setMethod("fitIndices", "regularizedSEMMixedPenalty", function(object) {
  
  fits <- object@fits
  
  usesLikelihood <- any(!is.na(fits$m2LL))
  
  if(usesLikelihood){
    multiGroup <- !is(object@inputArguments$lavaanModel, "lavaan")
    
    if(!multiGroup){
      dataset <- lavInspect(object@inputArguments$lavaanModel, "data")
      # remove empty rows:
      if(any(apply(dataset,1,function(x) all(is.na(x))))){
        warning("Your data set has rows where all observations are missing. lessSEM will",
                "remove those rows, but it is recommended to do so before fitting the models.")
        dataset <- dataset[!apply(dataset,1,function(x) all(is.na(x))),,drop = FALSE]
      }
      sampstats <- lavInspect(object@inputArguments$lavaanModel, "sampstat")
      N <- nrow(dataset)
    }
    
    # fit indices
    fits$AIC <- AIC(object)$AIC
    fits$BIC <- BIC(object)$BIC
    
    # The following variants of the AIC are adapted from here:
    # https://search.r-project.org/CRAN/refmans/AICcmodavg/html/AICc.html
    if(!multiGroup){
      fits$AICc <- fits$m2LL + 2*fits$nonZeroParameters * (N/(N - fits$nonZeroParameters - 1))
      
      # Chi^2
      
      if(is.null(sampstats$mean)){
        sampstats$mean <- apply(dataset, 2, mean, na.rm = TRUE)
        satPar <- nrow(sampstats$cov)*(ncol(sampstats$cov)+1)/2
      }else{
        satPar <- nrow(sampstats$cov)*(ncol(sampstats$cov)+1)/2 + length(sampstats$mean)
      }
      
      saturatedFit <- -2*sum(apply(dataset, 1, function(x) mvtnorm::dmvnorm(x = x[!is.na(x), drop = FALSE], 
                                                                            mean = sampstats$mean[!is.na(x), drop = FALSE], 
                                                                            sigma = sampstats$cov[!is.na(x), !is.na(x), drop = FALSE], 
                                                                            log = TRUE))
      )
      
      fits$chisq <- fits$m2LL - saturatedFit
      fits$df <- satPar - fits$nonZeroParameters
      
      # RMSEA
      # degrees of freedom
      lambda <- fits$chisq - fits$df
      # Note: lavaan uses df*N instead of df*(N-1)!
      N <- nrow(dataset)
      
      fits$rmsea <- 0
      fits$rmsea[lambda >= 0] <- sqrt(lambda[lambda>=0] / (fits$df[lambda>=0] * N))
    }
  }
  
  return(fits)
  
})

#' estimates
#' 
#' @param object object of class regularizedSEMMixedPenalty
#' @param criterion fit index (e.g., AIC) used to select the parameters
#' @param transformations boolean: Should transformations be returned?
#' @return returns a matrix with estimates
#' @export
setMethod("estimates", "regularizedSEMMixedPenalty", function(object, criterion = NULL, transformations = FALSE) {
  
  if(transformations)
    return(cbind(
      coef(object, criterion = criterion)@estimates,
      coef(object, criterion = criterion)@transformations)
    )
  
  return(coef(object, criterion = criterion)@estimates)
  
})