#### This script provides functions to perform approximate cross-validation for different types of models ####

#' aCV4regularizedSEM
#' 
#' approximate cross-validation for models of class regularizedSEM. These models can be fit with regularizedSEM() (see ?regularizedSEM)
#' in this package.
#' 
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param eps controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.
#' @export
aCV4regularizedSEM <- function(regularizedSEM, k, eps = 1e-4){
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  
  data <- try(lavaan::lavInspect(regularizedSEM$inputArguments$lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  aCVSEM <- SEMFromLavaan(lavaanModel = regularizedSEM$inputArguments$lavaanModel, rawData = data, transformVariances = TRUE)
  
  # extract tuning parameters
  penalty <- regularizedSEM$inputArguments$penalty
  lambdas <- regularizedSEM$inputArguments$lambdas
  alphas <- regularizedSEM$inputArguments$alphas
  
  if(penalty == "elasticNet"){
    tuningParameters <- expand.grid(lambda = lambdas, alpha = alphas)
  }else{
    tuningParameters <- expand.grid(lambda = lambdas)
  }
  
  if(penalty == "elasticNet") {
    coln <- paste0("lambda=",tuningParameters$lambda, ";alpha=", tuningParameters$alpha)
  }else{
    coln <- paste0("lambda=",tuningParameters$lambda)
  }
  
  regularizedParameters <- regularizedSEM$parameters
  regularizedParameterLabels <- regularizedSEM$inputArguments$regularizedParameterLabels
  
  aCVs <- matrix(NA,
                 nrow = k, 
                 ncol = nrow(regularizedParameters),
                 dimnames = list(
                   paste0("sample", 1:k),
                   coln
                 ))
  
  progressbar = txtProgressBar(min = 0,
                               max = nrow(regularizedParameters), 
                               initial = 0, 
                               style = 3)
  
  for(ro in 1:nrow(tuningParameters)){
    
    if(is.null(tuningParameters$alpha)){
      
      lambda <- tuningParameters$lambda[ro]
      selectPars <- paste0("lambda=",lambda)
      
    }else{
      
      lambda <- tuningParameters$alpha[ro]
      alpha <- tuningParameters$alpha[ro]
      selectPars <- paste0("lambda=",lambda, ";alpha=", alpha)
      
    }
    
    aCVSEM <- setParameters(SEM = aCVSEM,
                            labels = colnames(regularizedParameters),
                            values = regularizedParameters[selectPars,],
                            raw = FALSE)
    aCVSEM <- fit(aCVSEM)
    
    if(penalty == "elasticNet"){
      penaltyFunctionArguments <- list(
        "regularizedParameterLabels" = regularizedParameterLabels,
        "lambda" = lambda,
        "alpha" = alpha,
        "eps" = eps
      )
      aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                  k = k,
                                  individualPenaltyFunction = smoothElasticNet, 
                                  individualPenaltyFunctionGradient = smoothElasticNetGradient,
                                  individualPenaltyFunctionHessian = smoothElasticNetHessian,
                                  raw = FALSE, 
                                  penaltyFunctionArguments = penaltyFunctionArguments)
      
    }
    if(penalty == "ridge"){
      penaltyFunctionArguments <- list(
        "regularizedParameterLabels" = regularizedParameterLabels,
        "lambda" = lambda,
      )
      aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                  k = k,
                                  individualPenaltyFunction = ridge, 
                                  individualPenaltyFunctionGradient = ridgeGradient,
                                  individualPenaltyFunctionHessian = ridgeHessian,
                                  raw = FALSE, 
                                  penaltyFunctionArguments = penaltyFunctionArguments)
      
    }
    if(penalty == "lasso"){
      penaltyFunctionArguments <- list(
        "regularizedParameterLabels" = regularizedParameterLabels,
        "lambda" = lambda,
        "eps" = eps
      )
      aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                  k = k,
                                  individualPenaltyFunction = smoothLASSO, 
                                  individualPenaltyFunctionGradient = smoothLASSOGradient,
                                  individualPenaltyFunctionHessian = smoothLASSOHessian,
                                  raw = FALSE, 
                                  penaltyFunctionArguments = penaltyFunctionArguments)
    }
    
    aCVs[paste0("sample", 1:k),selectPars] <- aCV$leaveOutFits
    
    setTxtProgressBar(progressbar,which(colnames(aCVs) == selectPars))
    
  }
  
  return(list(
    "lambda" = lambdas,
    "alpha" = alphas,
    "parameters" = regularizedParameters,
    "leaveOutFits" = aCVs
  )
  )
}

#' aCV4lavaan
#' 
#' approximate cross-validation for models of class lavaan. 
#' 
#' @param lavaanModel model of class lavaan
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param raw controls if the cross-validation should use the internal transformations of aCV4SEM. aCV4SEM will use an exponential function for all variances to 
#' avoid negative variances. This can result in better sub-group parameters, but may not be necessary and will also result in more difficult to interpret parameters.
#' @export
aCV4lavaan <- function(lavaanModel, 
                       k, 
                       raw = FALSE){
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  data <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  aCVSEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = data, transformVariances = TRUE)
  
  return(aCV4SEM:::smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                        k = k,
                                        individualPenaltyFunction = NULL, 
                                        individualPenaltyFunctionGradient = NULL,
                                        individualPenaltyFunctionHessian = NULL,
                                        raw = raw, 
                                        penaltyFunctionArguments = NULL))
}

#' aCV4regsem
#' 
#' approximate cross-validation for models of class regsem. 
#' 
#' @param regsemModel model of class regsem
#' @param lavaanModel model of class lavaan. This must be the same model used to set up the regsem model!
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param penalty which penalty was used in regsem? Currently available are: penalty = "lasso", penalty = "ridge", or penalty = "elasticNet"
#' @param lambda value of the tuning parameter lambda. The tuning parameter alpha for the elastic net will be exported from the regsemModel
#' @param eps controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.
#' @export
aCV4regsem <- function(regsemModel, lavaanModel, k, penalty, lambda, eps = 1e-4){
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  if(!penalty %in% c("lasso", "ridge", "elasticNet")) stop("approximateCrossValidationCvregsem currently only supports lasso, ridge, or elastic net as penalty functions")
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  data <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  if(!is(regsemModel, "regsem")) stop("regsemModel must be of class regsem.")
  
  aCVSEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = data, transformVariances = TRUE)
  
  regsemParameters <- regsem2LavaanParameters(regsemModel = regsemModel, lavaanModel = lavaanModel)
  regularizedParameterLabels <- names(regsemParameters)[regsemModel$pars_pen]
  
  message("aCV4SEM assumes that the following parameters were regularized: ", 
          paste0(regularizedParameterLabels, collapse = ", "), 
          ". Please make sure that this is correct.")
  
  # extract tuning parameters
  if(penalty == "elasticNet") {
    if(is.null(regsemModel$call$alpha)) stop("Cannot find the tuning parameter alpha in regsemModel. Please make sure to call cv_regsem with alpha explicitly specified")
    alpha <- regsemModel$call$alpha # cv_regsem cannot use a grid of alphas
  }
  
  aCVSEM <- setParameters(SEM = aCVSEM,
                          labels = names(regsemParameters),
                          values = regsemParameters,
                          raw = FALSE)
  aCVSEM <- fit(aCVSEM)
  
  if(penalty == "elasticNet"){
    penaltyFunctionArguments <- list(
      "regularizedParameterLabels" = regularizedParameterLabels,
      "lambda" = lambda,
      "alpha" = alpha,
      "eps" = eps
    )
    aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                k = k,
                                individualPenaltyFunction = smoothElasticNet, 
                                individualPenaltyFunctionGradient = smoothElasticNetGradient,
                                individualPenaltyFunctionHessian = smoothElasticNetHessian,
                                raw = FALSE, 
                                penaltyFunctionArguments = penaltyFunctionArguments)
    
  }
  if(penalty == "ridge"){
    penaltyFunctionArguments <- list(
      "regularizedParameterLabels" = regularizedParameterLabels,
      "lambda" = lambda,
    )
    aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                k = k,
                                individualPenaltyFunction = ridge, 
                                individualPenaltyFunctionGradient = ridgeGradient,
                                individualPenaltyFunctionHessian = ridgeHessian,
                                raw = FALSE, 
                                penaltyFunctionArguments = penaltyFunctionArguments)
    
  }
  if(penalty == "lasso"){
    penaltyFunctionArguments <- list(
      "regularizedParameterLabels" = regularizedParameterLabels,
      "lambda" = lambda,
      "eps" = eps
    )
    aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                k = k,
                                individualPenaltyFunction = smoothLASSO, 
                                individualPenaltyFunctionGradient = smoothLASSOGradient,
                                individualPenaltyFunctionHessian = smoothLASSOHessian,
                                raw = FALSE, 
                                penaltyFunctionArguments = penaltyFunctionArguments)
  }
  
  
  
  
  return(list(
    "lambda" = lambda,
    "parameters" = regsemParameters,
    "leaveOutFits" = aCV$leaveOutFits,
    "aCVObject" = aCV)
  )
}

#' aCV4cv_regsem
#' 
#' approximate cross-validation for models of class cv_regsem. 
#' 
#' @param cvregsemModel model of class cv_regsem
#' @param lavaanModel model of class lavaan. This must be the same model used to set up the regsem model!
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param penalty which penalty was used in regsem? Currently available are: penalty = "lasso", penalty = "ridge", or penalty = "elasticNet"
#' @param exact boolean: If set to TRUE, an exact optimization based on GLMNET will be used to compute the sub-group parameters.
#' @param eps Only relevant if exact = FALSE. Controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.
#' @export
aCV4cv_regsem <- function(cvregsemModel, lavaanModel, k, penalty, exact = FALSE, eps = 1e-4){
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  if(!penalty %in% c("lasso", "ridge", "elasticNet")) stop("approximateCrossValidationCvregsem currently only supports lasso, ridge, or elastic net as penalty functions")
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  data <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  if(!is(cvregsemModel, "cvregsem")) stop("cvregsemModel must be of class regsem.")
  
  aCVSEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = data, transformVariances = TRUE)
  
  regsemParameters <- cvregsem2LavaanParameters(cvregsemModel = cvregsemModel, lavaanModel = lavaanModel)
  regularizedParameterLabels <- colnames(regsemParameters)[cvregsemModel$pars_pen]
  
  message("aCV4SEM assumes that the following parameters were regularized: ", 
          paste0(regularizedParameterLabels, collapse = ", "), 
          ". Please make sure that this is correct.")
  
  # extract tuning parameters
  lambdas <- cvregsemModel$fits[,"lambda"]
  if(penalty == "elasticNet") {
    if(is.null(cvregsemModel$call$alpha)) stop("Cannot find the tuning parameter alpha in cvregsemModel. Please make sure to call cv_regsem with alpha explicitly specified")
    alphas <- rep(cvregsemModel$call$alpha, length(lambdas)) # cv_regsem cannot use a grid of alphas
  }
  
  if(penalty == "elasticNet") {
    coln <- paste0("lambda=",lambdas, "; alpha=", alphas)
  }else{
    coln <- paste0("lambda=",lambdas)
  }
  aCVs <- matrix(NA,
                 nrow = k, 
                 ncol = nrow(regsemParameters),
                 dimnames = list(
                   paste0("sample", 1:k),
                   coln
                 ))
  
  progressbar = txtProgressBar(min = 0,
                               max = nrow(regsemParameters), 
                               initial = 0, 
                               style = 3)
  
  for(p in 1:nrow(regsemParameters)){
    if(anyNA(regsemParameters[p,])) next
    aCVSEM <- setParameters(SEM = aCVSEM,
                            labels = colnames(regsemParameters),
                            values = regsemParameters[p,],
                            raw = FALSE)
    aCVSEM <- fit(aCVSEM)
    
    if(penalty == "elasticNet"){
      penaltyFunctionArguments <- list(
        "regularizedParameterLabels" = regularizedParameterLabels,
        "lambda" = lambdas[p],
        "alpha" = alphas[p],
        "eps" = eps
      )
      aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                  k = k,
                                  individualPenaltyFunction = smoothElasticNet, 
                                  individualPenaltyFunctionGradient = smoothElasticNetGradient,
                                  individualPenaltyFunctionHessian = smoothElasticNetHessian,
                                  raw = FALSE, 
                                  penaltyFunctionArguments = penaltyFunctionArguments)
      
    }
    if(penalty == "ridge"){
      penaltyFunctionArguments <- list(
        "regularizedParameterLabels" = regularizedParameterLabels,
        "lambda" = lambdas[p],
      )
      aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                  k = k,
                                  individualPenaltyFunction = ridge, 
                                  individualPenaltyFunctionGradient = ridgeGradient,
                                  individualPenaltyFunctionHessian = ridgeHessian,
                                  raw = FALSE, 
                                  penaltyFunctionArguments = penaltyFunctionArguments)
      
    }
    if(penalty == "lasso"){
      
      if(exact){
        penaltyFunctionArguments <- list(
          "regularizedParameterLabels" = regularizedParameterLabels,
          "lambda" = lambdas[p]
        )
        aCV <- GLMNETACVRcpp_SEMCpp(SEM = aCVSEM,
                                    k = k,
                                    penalty = penalty,
                                    raw = TRUE,
                                    penaltyFunctionArguments = penaltyFunctionArguments)
      }else{
        penaltyFunctionArguments <- list(
          "regularizedParameterLabels" = regularizedParameterLabels,
          "lambda" = lambdas[p],
          "eps" = eps
        )
        aCV <- smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                    k = k,
                                    individualPenaltyFunction = smoothLASSO, 
                                    individualPenaltyFunctionGradient = smoothLASSOGradient,
                                    individualPenaltyFunctionHessian = smoothLASSOHessian,
                                    raw = FALSE, 
                                    penaltyFunctionArguments = penaltyFunctionArguments)
        
      }
      
    }
    
    aCVs[paste0("sample", 1:k),p] <- aCV$leaveOutFits
    
    setTxtProgressBar(progressbar,p)
  }
  
  return(list(
    "lambda" = lambdas,
    "regsemParameters" = regsemParameters,
    "leaveOutFits" = aCVs
  )
  )
}
