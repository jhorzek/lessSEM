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

optimizeSEMLASSO <- function(SEM, alternativeStartingValues = NULL, regularizedParameters, lambda, raw = TRUE, ...){
  parameterBackup <- getParameters(SEM, raw = TRUE) # save in case something brakes
  
  # get parameters
  parameters <- getParameters(SEM, raw = raw)
  
  # check regularizedParameters
  if(any(!regularizedParameters %in% names(parameters))) stop(paste0("Could not find parameter ", 
                                                                     paste0(regularizedParameters[!regularizedParameters %in% names(parameters)], collapse = ", "))
  )
  
  # check if variances are regularized and the parameters are raw
  if(raw){
    parameterTable <- SEM$getParameters()
    for(regularizedParmeter in regularizedParameters){
      if(parameterTable$location[parameterTable$label == regularizedParmeter] != "Smatrix") next
      if(parameterTable$row[parameterTable$label == regularizedParmeter] == parameterTable$col[parameterTable$label == regularizedParmeter]){
        warning("Regularizing raw_value of variances (raw = TRUE). The actual variances are given by var = exp(raw_value).")
      }
    }
  }
  
  Nobs <- nrow(SEM$rawData)
  
  penaltyFunction <- function(par, SEM, raw, lambda, regularizedParameters, Nobs){
    return(Nobs*lambda*sum(sqrt((par[regularizedParameters])^2 + 1e-10)))
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

    gradientsM2LL <- derivativeFunction(par, SEM, raw)
    if(any(class(gradientsM2LL) == "try-error")) return(rep(9999, length(par)))
    # add derivative of regularization function
    
    gradientsM2LL[regularizedParameters] <- gradientsM2LL[regularizedParameters] + 
      Nobs*lambda*par[regularizedParameters]/(sqrt(par[regularizedParameters]^2 + 1e-10))
    
    return(gradientsM2LL)
  }
  
  optimized <- optim(par = parameters, 
                     fn = fitfun, 
                     gr = computeRegularizedDerivativesF, 
                     penaltyFunction = penaltyFunction,
                     SEM = SEM, 
                     raw = raw,
                     lambda = lambda, 
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
                       lambda = lambda, 
                       regularizedParameters = regularizedParameters, 
                       Nobs = Nobs, 
                       method = "BFGS")
    
    
    
    SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
    SEM <- try(fit(SEM), silent = TRUE)
  }
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}