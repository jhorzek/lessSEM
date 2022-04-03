optimizeSEM <- function(SEM, raw = TRUE, ...){
  # get parameters and check for variances -> these will get a lower bound
  parameters <- getParameters(SEM, raw = raw)
  
  fitfun <- function(par, 
                     SEM,
                     raw){
    SEM <- setParameters(SEM, names(par), par, raw = raw)
    SEM <- try(fit(SEM), silent = TRUE)
    
    if(any(class(SEM) == "try-error") || is.na(SEM$m2LL)) return(999999999999)
    
    return(SEM$m2LL)
    
  }
  
  fitDerivative <- function(par, SEM, raw){
    tryDerivs <- try(computeAnalyticGradients(par, SEM, raw = raw), silent = TRUE)
    if(any(class(tryDerivs) == "try-error")) return(rep(9999, length(par)))
    return(tryDerivs)
  }
  
  optimized <- optim(par = parameters, 
                     fn = fitfun, 
                     gr = fitDerivative, 
                     SEM = SEM,
                     raw = raw,
                     method = "BFGS",
                     ...)
  SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(fit(SEM), silent = TRUE)
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}

optimizeLassoSEM <- function(SEM, alternativeStartingValues = NULL, regularizedParameters, lambda, raw = TRUE, ...){
  SEM_ <- SEM # save in case something brakes
  # get parameters and check for variances -> these will get a lower bound
  parameters <- getParameters(SEM, raw = raw)
  parameterTable <- SEM$model$parameters
  variances <- parameterTable$label[parameterTable$row == parameterTable$col & parameterTable$location == "Smatrix"]
  
  # check regularizedParameters
  if(any(!regularizedParameters %in% names(parameters))) stop(paste0("Could not find parameter ", 
                                                                     regularizedParameters[!regularizedParameters %in% names(parameters)]))
  
  Nobs <- nrow(SEM$data$rawData)
  
  penaltyFunction <- function(par, SEM, raw, lambda, regularizedParameters, Nobs){
    SEM <- setParameters(SEM, names(par), par, raw = raw) # in case the parameters were raw
    parTransformed <- getParameters(SEM, raw = FALSE)
    return(lambda*sum(sqrt((parTransformed[regularizedParameters])^2 + 1e-10)))
  }
  
  fitfun <- function(par, 
                     SEM, 
                     raw,
                     penaltyFunction,
                     lambda, 
                     regularizedParameters, 
                     Nobs){
    SEM <- setParameters(SEM, names(par), par, raw = raw)
    SEM <- try(fit(SEM), silent = TRUE)
    
    if(any(class(SEM) == "try-error") || is.na(SEM$m2LL)) {
      return(9999999999)
      }
    
    ff <- SEM$m2LL + penaltyFunction(par = par, 
                                     SEM = SEM, 
                                     raw = raw, 
                                     lambda = lambda, 
                                     regularizedParameters = regularizedParameters, 
                                     Nobs = Nobs)
    
    return(ff)
    
  }
  
  computeRegularizedDerivativesF <- function(par, 
                                             SEM, 
                                             raw,
                                             penaltyFunction,
                                             penaltyFunctionDerivative = NULL, 
                                             lambda, 
                                             regularizedParameters,
                                             Nobs){
    SEM <- setParameters(SEM, names(par), par, raw = raw) # in case the parameters were raw
    parTransformed <- getParameters(SEM, raw = FALSE)
    gradientsM2LL <- try(computeAnalyticGradients(par = par, SEM = SEM, raw = raw))
    if(any(class(gradientsM2LL) == "try-error")) return(rep(9999, length(par)))
    # add derivative of regularization function
    
    gradientsM2LL[regularizedParameters] <- gradientsM2LL[regularizedParameters] + 
      lambda*parTransformed[regularizedParameters]/(sqrt(parTransformed[regularizedParameters]^2 + 1e-10))
    
    return(gradientsM2LL)
  }
  
  optimized <- optim(par = parameters, 
                     fn = fitfun, 
                     gr = computeRegularizedDerivativesF, 
                     penaltyFunction = penaltyFunction,
                     SEM = SEM, 
                     raw = raw,
                     lambda = Nobs*lambda, 
                     regularizedParameters = regularizedParameters, 
                     Nobs = Nobs, 
                     method = "BFGS")
  


  SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(fit(SEM), silent = TRUE)
  
  if(any(class(SEM) == "try-error" || is.na(SEM$m2LL)) && !is.null(alternativeStartingValues)){
    print("Trying alternative starting values")
    optimized <- optim(par = alternativeStartingValues, 
                       fn = fitfun, 
                       gr = computeRegularizedDerivativesF, 
                       penaltyFunction = penaltyFunction,
                       SEM = SEM_, 
                       raw = raw,
                       lambda = Nobs*lambda, 
                       regularizedParameters = regularizedParameters, 
                       Nobs = Nobs, 
                       method = "BFGS")
    
    
    
    SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
    SEM <- try(fit(SEM), silent = TRUE)
  }
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}