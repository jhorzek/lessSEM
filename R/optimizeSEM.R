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

optimizeRegularizedSEM <- function(SEM, raw = TRUE, individualPenaltyFunction, controlOptim = NULL, method = "BFGS", ...){
  
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

optimizeSEMLASSO <- function(SEM, regularizedParameters, lambda, eps = 1e-4, raw = TRUE, ...){

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
  
  individualPenaltyFunction <- function(par, SEM, raw, lambda, regularizedParameters, eps){
    return(lambda*sum(sqrt((par[regularizedParameters])^2 + eps)))
  }
  
  out <- optimizeRegularizedSEM(SEM = SEM, 
                                raw = raw, 
                                individualPenaltyFunction = individualPenaltyFunction, 
                                lambda = lambda,
                                regularizedParameters = regularizedParameters, 
                                eps = eps,
                                method = "BFGS")
  
  return(out)
}