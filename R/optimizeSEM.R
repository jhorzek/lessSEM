optimizeSEM <- function(SEM, raw = TRUE, ...){
  
  parameters <- getParameters(SEM, raw = raw)
  
  optimized <- optim(par = parameters, 
                     fn = fitFunction, 
                     gr = derivativeFunction, 
                     SEM = SEM,
                     raw = raw,
                     method = "BFGS",
                     ...)
  SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(fit(SEM), silent = TRUE)
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}

optimizeCustomRegularizedSEM <- function(SEM, individualPenaltyFunction, controlOptim = NULL, method = "BFGS", raw = TRUE, ...){
  
  additionalArguments <- list(...)
  
  parameters <- getParameters(SEM, raw = raw)
  
  N <- nrow(SEM$rawData)
  
  fitFun <- function(par, SEM, raw, individualPenaltyFunction, ...){
    m2LL <- fitFunction(par = par, SEM = SEM, raw = raw)
    if(m2LL == 9999999999999) return(m2LL)
    
    m2LLRegularized <- m2LL + N*individualPenaltyFunction(par, ...)
    
    return(m2LLRegularized)
  }
  
  gradFun <- function(par, SEM, raw, individualPenaltyFunction, ...){
    gradients <- derivativeFunction(par = par, SEM = SEM, raw = raw)
    if(all(gradients == 9999999999999)) return(gradients)
    
    gradients <- gradients + N*numDeriv::grad(func = individualPenaltyFunction, 
                                              x = par, 
                                              ...)
    
    return(gradients)
  }
  
  if(is.null(controlOptim)) {
    optimized <- optim(par = parameters, 
                       fn = fitFun, 
                       gr = gradFun, 
                       SEM = SEM,
                       raw = raw,
                       individualPenaltyFunction = individualPenaltyFunction,
                       method = method,
                       ...)
  }else{
    optimized <- optim(par = parameters, 
                       fn = fitFun, 
                       gr = gradFun, 
                       SEM = SEM,
                       raw = raw,
                       individualPenaltyFunction = individualPenaltyFunction,
                       method = method,
                       control = controlOptim,
                       ...)
  }
  
  SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(fit(SEM), silent = TRUE)
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}

optimizeRegularizedSEM <- function(lavaanModel, regularizedParameterLabels, penalty, lambda, alpha = NULL, eps = 1e-4, controlOptim = NULL, method = "BFGS", raw = TRUE){
  inputArguments <- as.list(environment())
  
  if(eps < 1e-5) warning("Setting eps too small may result in very ragged parameter estimates. We recommend trying different values of eps.")
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
  
  if(penalty == "ridge") {
    result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                           raw = raw, 
                                           individualPenaltyFunction = ridge, 
                                           lambda = lambda,
                                           regularizedParameterLabels = regularizedParameterLabels, 
                                           controlOptim = controlOptim,
                                           method = method)
    return(
      list(
        "parameters" = getParameters(result$SEM, raw = FALSE),
        "result" = result,
        "inputArguments" = inputArguments)
    )
  }
  
  if(penalty == "lasso"){ 
    result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                           raw = raw, 
                                           individualPenaltyFunction = smoothLASSO, 
                                           lambda = lambda,
                                           regularizedParameterLabels = regularizedParameterLabels, 
                                           eps = eps,
                                           controlOptim = controlOptim,
                                           method = method)
    return(
      list(
        "parameters" = getParameters(result$SEM, raw = FALSE),
        "result" = result,
        "inputArguments" = inputArguments)
    )
  }
  
  if(penalty == "adaptiveLasso") {
    if(length(lambda) != length(regularizedParameterLabels)) stop("adaptiveLasso requires a lambda value for each regularized parameter.")
    if(any(!names(lambda) %in% regularizedParameterLabels)){stop("names of lambda do not match the regularizedParameterLabels vector.")}
    result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                           raw = raw, 
                                           individualPenaltyFunction = smoothAdaptiveLASSO, 
                                           lambda = lambda,
                                           regularizedParameterLabels = regularizedParameterLabels, 
                                           controlOptim = controlOptim,
                                           method = method)
    return(
      list(
        "parameters" = getParameters(result$SEM, raw = FALSE),
        "result" = result,
        "inputArguments" = inputArguments)
    )
  }
  
  if(penalty == "elasticNet") {
    if(is.null(alpha)) stop("elasticNet requires specification of alpha.")
    result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                           raw = raw, 
                                           individualPenaltyFunction = smoothElasticNet, 
                                           lambda = lambda,
                                           alpha = alpha,
                                           regularizedParameterLabels = regularizedParameterLabels, 
                                           eps = eps,
                                           controlOptim = controlOptim,
                                           method = method)
    return(
      list(
        "parameters" = getParameters(result$SEM, raw = FALSE),
        "result" = result,
        "inputArguments" = inputArguments)
    )
  }
  stop("Something went terribly wrong...")
}

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