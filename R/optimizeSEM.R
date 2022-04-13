#' optimizeSEM
#' 
#' optimize internal SEM representation with optim() function
#' 
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of aCV4SEM is used.
#' @param method optimizer method. See ?optim
#' @param control control optimizer. See ?optim
#' @param hessian Should the optimizer return the hessian?
#' @export
optimizeSEM <- function(SEM, raw = TRUE, method = "BFGS", control = list(), hessian = FALSE){
  
  parameters <- getParameters(SEM, raw = raw)
  
  optimized <- optim(par = parameters, 
                     fn = fitFunction, 
                     gr = derivativeFunction, 
                     SEM = SEM,
                     raw = raw,
                     method = method,
                     control = control,
                     hessian = hessian)
  SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(fit(SEM), silent = TRUE)
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}

#' optimizeCustomRegularizedSEM
#' 
#' optimize a SEM with custom penalty function.
#' 
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument and the penaltyFunctionArguments as second argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionGradient gradients of the penalty function. Function should take the current parameter values as first argument and the penaltyFunctionArguments as second argument and 
#' return a vector of the same length as parameters. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionHessian Hessian of the penalty function. Function should take the current parameter values as first argument and the penaltyFunctionArguments as second argument and 
#' return a matrix with (length as parameters)^2 number of elements. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param raw controls if the internal transformations of aCV4SEM is used.
#' @param method optimizer method. See ?optim
#' @param control control optimizer. See ?optim
#' @param hessian Should the optimizer return the hessian?
#' @export
optimizeCustomRegularizedSEM <- function(SEM, 
                                         individualPenaltyFunction, 
                                         individualPenaltyFunctionGradient, 
                                         penaltyFunctionArguments,
                                         raw = TRUE, 
                                         method = "BFGS", 
                                         control = list(), 
                                         hessian = FALSE){
  
  
  parameters <- getParameters(SEM, raw = raw)
  
  N <- nrow(SEM$rawData)
  
  fitFun <- function(par, SEM, raw, individualPenaltyFunction, individualPenaltyFunctionGradient, penaltyFunctionArguments){
    m2LL <- fitFunction(par = par, SEM = SEM, raw = raw)
    if(m2LL == 9999999999999) return(m2LL)
    m2LLRegularized <- m2LL 
    for(i in 1:N){
      m2LLRegularized <- m2LLRegularized + individualPenaltyFunction(par, penaltyFunctionArguments)
    }
    return(m2LLRegularized)
  }
  
  gradFun <- function(par, SEM, raw, individualPenaltyFunction, individualPenaltyFunctionGradient, penaltyFunctionArguments){
    gradients <- derivativeFunction(par = par, SEM = SEM, raw = raw)
    if(all(gradients == 9999999999999)) return(gradients)
    
    for(i in 1:N){
      gradients <- gradients + individualPenaltyFunctionGradient(par, penaltyFunctionArguments)
    }
    
    return(gradients)
  }
  
  
  optimized <- optim(par = parameters, 
                     fn = fitFun, 
                     gr = gradFun, 
                     SEM = SEM,
                     raw = raw,
                     individualPenaltyFunction = individualPenaltyFunction,
                     individualPenaltyFunctionGradient = individualPenaltyFunctionGradient,
                     penaltyFunctionArguments = penaltyFunctionArguments,
                     method = method,
                     control = control,
                     hessian = hessian)
  
  SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(fit(SEM), silent = TRUE)
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}

#' optimizeRegularizedSEM
#' 
#' optimize a SEM with regularized penalty function. Uses smooth approximations for non-differentiable penalty functions.
#' 
#' @param lavaanModel model of class lavaan 
#' @param regularizedParameterLabels labels of regularized parameters
#' @param penalty which penalty should be used? Available are "ridge", "lasso", "adaptiveLasso", and "elasticNet"
#' @param lambdas vector with lambda values. Higher values = higher penalty
#' @param alphas 0<alpha<1. only required for elastic net. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge
#' @param adaptiveLassoWeights vector with weights for adaptive LASSO. Set to NULL if not using adaptive LASSO
#' @param eps controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.
#' @param raw controls if the internal transformations of aCV4SEM is used.
#' @param method optimizer method. See ?optim
#' @param control control optimizer. See ?optim
#' @param hessian Should the optimizer return the hessian?
#' @export
optimizeRegularizedSEM <- function(lavaanModel, 
                                   regularizedParameterLabels,
                                   penalty, 
                                   lambdas, 
                                   alphas = NULL, 
                                   adaptiveLassoWeights = NULL,
                                   raw = TRUE,
                                   method = "BFGS", 
                                   control = list(), 
                                   hessian = FALSE
){
  
  inputArguments <- as.list(environment())
  
  if(!penalty %in% c("lasso", "ridge", "adaptiveLasso", "elasticNet")) stop("Currently supported penalty functions are: lasso, ridge, adaptiveLasso, and elasticNet")
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  
  SEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = rawData, transformVariances = TRUE)
  
  # get parameters
  parameters <- getParameters(SEM, raw = raw)
  checkRegularizedParameters(parameters = parameters, 
                             regularizedParameterLabels = regularizedParameterLabels,
                             SEM$getParameters(), 
                             raw = raw)
  
  initialHessian <- getHessian(SEM = SEM, raw = TRUE)
  initialHessian4Optimizer <- initialHessian
  
  if(penalty == "adaptiveLasso" && is.null(adaptiveLassoWeights)){
    adaptiveLassoWeights <- abs(parameters)^(-1)
  }else{
    adaptiveLassoWeights <- rep(1, length(parameters))
    names(adaptiveLassoWeights) <- names(parameters)
  }
  inputArguments$adaptiveLassoWeights <- adaptiveLassoWeights
  
  if(!is.null(alphas) && penalty != "elasticNet") {stop("non-null alpha parameter only valid for elasticNet")}
  if(is.null(alphas)) {
    if(penalty %in% c("lasso", "adaptiveLasso")) {
      alphas <- 1
    }else if(penalty == "ridge"){
      alphas <- 0
    }
  }
  
  parameterEstimates <- matrix(NA, 
                               nrow = length(lambdas)*length(alphas), 
                               ncol = length(parameters))
  colnames(parameterEstimates) <- names(parameters)
  
  rownames(parameterEstimates) <- paste0("lambda=", lambdas, ";", rep(paste0("alpha=", alphas), each = length(lambdas)))
  
  Hessians <- array(NA, 
                    dim = c(length(parameters), length(parameters), nrow(parameterEstimates)),
                    dimnames = list(names(parameters),
                                    names(parameters),
                                    rownames(parameterEstimates))
  )
  it <- 0
  
  progressbar = txtProgressBar(min = 0, 
                               max = nrow(parameterEstimates), 
                               initial = it, 
                               style = 3)
  
  for(lambda in lambdas){
    
    for(alpha in alphas){
      
      whereTo <- paste0("lambda=",lambda, ";", "alpha=",alpha)
      
      result <- GLMNET(SEM = SEM, 
                       regularizedParameterLabels = regularizedParameterLabels, 
                       penalty = penalty, 
                       lambda = ifelse(penalty == "ridge", 2*lambda, lambda), # we are using the 
                       # elastic net implementation which takes .5*lambda*(1-alpha) with alpha = 0
                       # Setting lambda to 2*lambda makes sure that we are getting the ridge 
                       # estimates for lambda, not -5*lambda.
                       alpha = alpha, 
                       adaptiveLassoWeights = adaptiveLassoWeights,
                       initialHessian = initialHessian4Optimizer)
      
      newParameters <- getParameters(result$SEM, raw = FALSE)
      parameterEstimates[whereTo,] <- newParameters
      Hessians[,,whereTo] <- result$Hessian
      
      # set initial values for next iteration
      SEM <- setParameters(SEM = SEM, labels = names(newParameters), value = newParameters, raw = FALSE)
      SEM <- try(fit(SEM))
      if(is(SEM, "try-Error")){
        # reset
        warning("Fit for ", whereTo, " resultet in Error!")
        SEM <- SEMFromLavaan(lavaanModel = lavaanModel, rawData = rawData, transformVariances = TRUE)
        initialHessian4Optimizer <- NULL
      }else{
        initialHessian4Optimizer <- result$Hessian
      }
      
      it <- it + 1
      setTxtProgressBar(progressbar,it)
    }
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  results <- list(
    "parameters" = parameterEstimates,
    "internalOptimization" = internalOptimization,
    "inputArguments" = inputArguments)
  class(results) <- "regularizedSEM"
  return(
    results
  )
  
}

#' checkRegularizedParameters
#' 
#' internal function to check if the regularized parameters can be found in the model
#' 
#' @param parameters labeled vector with parameter values
#' @param regularizedParameterLabels labels of regularized parameters
#' @param parameterTable table with all parameters
#' @param raw controls if the internal transformations of aCV4SEM is used.
#' @export
checkRegularizedParameters <- function(parameters, regularizedParameterLabels, parameterTable, raw){
  
  
  # check regularizedParameterLabels
  if(any(!regularizedParameterLabels %in% names(parameters))) stop(paste0("Could not find parameter ", 
                                                                          paste0(regularizedParameterLabels[!regularizedParameterLabels %in% names(parameters)], collapse = ", "))
  )
  
  # check if variances are regularized and the parameters are raw
  if(raw){
    for(regularizedParmeter in regularizedParameterLabels){
      if(parameterTable$location[parameterTable$label == regularizedParmeter] != "Smatrix") next
      if(parameterTable$row[parameterTable$label == regularizedParmeter] == parameterTable$col[parameterTable$label == regularizedParmeter]){
        warning("Regularizing raw_value of variances (raw = TRUE). The actual variances are given by var = exp(raw_value).")
      }
    }
  }
}