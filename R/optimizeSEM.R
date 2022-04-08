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

optimizeCustomRegularizedSEM <- function(SEM, individualPenaltyFunction, individualPenaltyFunctionGradient, penaltyFunctionArguments, controlOptim = NULL, method = "BFGS", raw = TRUE, ...){
  
  additionalArguments <- list(...)
  
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
  
  if(is.null(controlOptim)) {
    optimized <- optim(par = parameters, 
                       fn = fitFun, 
                       gr = gradFun, 
                       SEM = SEM,
                       raw = raw,
                       individualPenaltyFunction = individualPenaltyFunction,
                       individualPenaltyFunctionGradient = individualPenaltyFunctionGradient,
                       penaltyFunctionArguments = penaltyFunctionArguments,
                       method = method
    )
  }else{
    optimized <- optim(par = parameters, 
                       fn = fitFun, 
                       gr = gradFun, 
                       SEM = SEM,
                       raw = raw,
                       individualPenaltyFunction = individualPenaltyFunction,
                       individualPenaltyFunctionGradient = individualPenaltyFunctionGradient,
                       penaltyFunctionArguments = penaltyFunctionArguments,
                       method = method,
                       control = controlOptim,
    )
  }
  
  SEM <- setParameters(SEM, names(optimized$par), optimized$par, raw = raw)
  SEM <- try(fit(SEM), silent = TRUE)
  
  return(list("SEM" = SEM,
              "optimizer" = optimized))
}

optimizeRegularizedSEM <- function(lavaanModel, 
                                   regularizedParameterLabels,
                                   penalty, 
                                   lambdas, 
                                   alphas = NULL, 
                                   eps = 1e-4, 
                                   controlOptim = NULL, 
                                   method = "BFGS", 
                                   raw = TRUE){
  
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
  
  if(penalty == "adaptiveLasso"){
    parameterWeights <- abs(parameters)^(-1)
  }
  
  if(!is.null(alphas) && penalty != "adaptiveLasso") {stop("non-null alpha parameter only valid for adaptiveLasso")}
  if(is.null(alphas)) alphas <- 0
  
  parameterEstimates <- matrix(NA, 
                              nrow = length(lambdas)*length(alphas), 
                              ncol = length(parameters))
  colnames(parameterEstimates) <- names(parameters)
  if(all(alphas == 0)){
    rownames(parameterEstimates) <- paste0("lambda=", lambdas)
  }else{
    rownames(parameterEstimates) <- paste0("lambda=", lambdas, ";", rep(paste0("alpha=", alphas), each = length(lambdas)))
  }
  
  optimResults <- vector("list", nrow(parameterEstimates))
  names(optimResults) <- rownames(parameterEstimates)
  
  progressbar = txtProgressBar(min = 0, 
                               max = nrow(parameterEstimates), 
                               initial = 0, 
                               style = 3)
  
  for(lambda in lambdas){
    
    for(alpha in alphas){
      
      if(penalty == "ridge") {
        result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                               raw = raw, 
                                               individualPenaltyFunction = ridge, 
                                               individualPenaltyFunctionGradient = ridgeGradient,
                                               penaltyFunctionArguments = list("lambda" = lambda,
                                                                               "regularizedParameterLabels" = regularizedParameterLabels
                                               ),
                                               controlOptim = controlOptim,
                                               method = method)
      }
      
      if(penalty == "lasso"){ 
        result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                               raw = raw, 
                                               individualPenaltyFunction = smoothLASSO, 
                                               individualPenaltyFunctionGradient = smoothLASSOGradient,
                                               penaltyFunctionArguments = list("lambda" = lambda,
                                                                               "regularizedParameterLabels" = regularizedParameterLabels,
                                                                               "eps" = eps
                                               ),
                                               controlOptim = controlOptim,
                                               method = method)
      }
      
      if(penalty == "adaptiveLasso") {
        if(length(lambda) != length(regularizedParameterLabels)) stop("adaptiveLasso requires a lambda value for each regularized parameter.")
        if(any(!names(lambda) %in% regularizedParameterLabels)){stop("names of lambda do not match the regularizedParameterLabels vector.")}
        result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                               raw = raw, 
                                               individualPenaltyFunction = smoothAdaptiveLASSO, 
                                               individualPenaltyFunctionGradient = smoothAdaptiveLASSOGradient,
                                               penaltyFunctionArguments = list("lambdas" = parameterWeights*lambda,
                                                                               "regularizedParameterLabels" = regularizedParameterLabels,
                                                                               "eps" = eps
                                               ),
                                               controlOptim = controlOptim,
                                               method = method)
      }
      
      if(penalty == "elasticNet") {
        if(is.null(alpha)) stop("elasticNet requires specification of alpha.")
        result <- optimizeCustomRegularizedSEM(SEM = SEM, 
                                               raw = raw, 
                                               individualPenaltyFunction = smoothElasticNet, 
                                               individualPenaltyFunctionGradient = smoothElasticNetGradient,
                                               penaltyFunctionArguments = list("lambda" = lambda,
                                                                               "alpha" = alpha,
                                                                               "regularizedParameterLabels" = regularizedParameterLabels,
                                                                               "eps" = eps
                                               ),
                                               controlOptim = controlOptim,
                                               method = method)
      }
      
      whereTo <- ifelse(all(alphas == 0), paste0("lambda=",lambda), paste0("lambda=",lambda, ";", "alpha=",alpha))
      optimResults[[whereTo]] <- result$optimized
      parameterEstimates[whereTo,] <- getParameters(result$SEM, raw = FALSE)
      
      setTxtProgressBar(progressbar,which(rownames(parameterEstimates) == whereTo))
    }
  }
  
  results <- list(
    "parameters" = parameterEstimates,
    "optimResults" = optimResults,
    "inputArguments" = inputArguments)
  class(results) <- "regularizedSEM"
  return(
    results
  )
  
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