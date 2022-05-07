#' regularizeSEMWithCustomPenalty
#' 
#' optimize a SEM with custom penalty functions.
#' 
#' @param lavaanModel model of class lavaan 
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionGradient gradients of the penalty function. Function should take the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' return a vector of the same length as parameters. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided. If set to NULL, numDeriv will be used.
#' @param individualPenaltyFunctionHessian Hessian of the penalty function. Function should take the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' return a matrix with (length as parameters)^2 number of elements. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided. If set to NULL, numDeriv will be used.
#' @param tuningParameters data.frame with tuning parameter values. Important: The function will iterate over the rows of these tuning parameters and pass them to your penalty function
#' @param penaltyFunctionArguments arguments passed to individualPenaltyFunction, individualPenaltyFunctionGradient, and individualPenaltyFunctionHessian
#' @param control option to set parameters of the optimizer
#' @export
regularizeSEMWithCustomPenalty <- function(lavaanModel, 
                                           individualPenaltyFunction,
                                           individualPenaltyFunctionGradient = NULL,
                                           individualPenaltyFunctionHessian = NULL,
                                           tuningParameters,
                                           penaltyFunctionArguments,
                                           control = aCV4SEM::controlQuasiNewtonBFGS()){
  
  if(is.null(individualPenaltyFunctionGradient)){
    message("Using numDeriv to approximate the gradient of the individualPenaltyFunction. This may result in slow optimization. We highly recommend that you pass an analytic solution.")
    individualPenaltyFunctionGradient <- aCV4SEM::genericGradientApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  if(is.null(individualPenaltyFunctionHessian)){
    message("Using numDeriv to approximate the Hessian of the individualPenaltyFunction. This may result in slow optimization. We highly recommend that you pass an analytic solution.")
    individualPenaltyFunctionHessian <- aCV4SEM::genericHessianApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  
  inputArguments <- as.list(environment())
  
  if(!is.data.frame(tuningParameters)) stop("tuningParameters must be a data frame (e.g., data.frame(lambda = seq(0,1,.1))).")
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  ### initialize model ####
  startingValues <- control$startingValues
  if(any(startingValues == "est")){
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "est")
  }else if(any(startingValues == "start")){
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start")
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(aCV4SEM::getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See aCV4SEM::getLavaanParameters(lavaanModel).")
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start", 
                                   fit = FALSE)
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(startingValues), value = startingValues, raw = FALSE)
    SEM <- try(aCV4SEM:::fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to regularizeSEMWithCustomPenalty. See ?controlQuasiNewtonBFGS for more information.")
  }
  
  # get parameters
  parameters <- aCV4SEM:::getParameters(SEM, raw = TRUE)
  
  initialHessian <- control$initialHessian
  if(is.matrix(initialHessian) && nrow(initialHessian) == length(parameters) && ncol(initialHessian) == length(parameters)){
    if(!all(rownames(initialHessian) %in% names(aCV4SEM::getLavaanParameters(lavaanModel))) ||
       !all(colnames(initialHessian) %in% names(aCV4SEM::getLavaanParameters(lavaanModel)))
    ) stop("initialHessian must have the parameter names as rownames and colnames. See aCV4SEM::getLavaanParameters(lavaanModel).")
  }else if(any(initialHessian == "compute")){
    initialHessian <- aCV4SEM:::getHessian(SEM = SEM, raw = TRUE)
  }else if(length(initialHessian) == 1 && is.numeric(initialHessian)){
    initialHessian <- diag(initialHessian,length(parameters))
    rownames(initialHessian) <- names(parameters)
    colnames(initialHessian) <- names(parameters)
  }else{
    stop("Invalid initialHessian passed to GLMNET. See ?controlGLMNET for more information.")
  }
  
  initialHessian4Optimizer <- initialHessian
  
  fits <- data.frame(
    "m2LL" = rep(NA,nrow(tuningParameters)),
    "regM2LL"= rep(NA,nrow(tuningParameters)),
    "convergence" = rep(NA,nrow(tuningParameters))
  )
  
  fits <- cbind(tuningParameters, fits)
  
  parameterEstimates <- as.data.frame(matrix(NA,nrow = nrow(tuningParameters), ncol = length(parameters)))
  colnames(parameterEstimates) <- names(parameters)
  parameterEstimates <- cbind(
    tuningParameters,
    parameterEstimates
  )
  
  if(control$saveHessian){
    Hessians <- list(
      "lambda" = tuningGrid$lambda,
      "alpha" = tuningGrid$alpha,
      "Hessian" = lapply(1:nrow(tuningGrid), 
                         matrix, 
                         data= NA, 
                         nrow=nrow(initialHessian), 
                         ncol=ncol(initialHessian))
    )
  }else{
    Hessians <- list(NULL)
  }
  
  if(control$verbose == 0){
    progressbar = txtProgressBar(min = 0, 
                                 max = nrow(tuningParameters), 
                                 initial = 0, 
                                 style = 3)
  }
  
  for(i in 1:nrow(tuningParameters)){
    
    if(control$verbose == 0){
      setTxtProgressBar(progressbar,i)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tuningParameters),"]\n"))
    }
    
    currentTuningParameters <- tuningParameters[i,,drop = FALSE]
    
    result <- aCV4SEM:::quasiNewtonBFGS(SEM = SEM, 
                                        individualPenaltyFunction = individualPenaltyFunction, 
                                        individualPenaltyFunctionGradient = individualPenaltyFunctionGradient, 
                                        individualPenaltyFunctionHessian = individualPenaltyFunctionHessian,
                                        currentTuningParameters = currentTuningParameters,
                                        penaltyFunctionArguments = penaltyFunctionArguments,
                                        initialHessian = initialHessian4Optimizer,
                                        stepSize = control$stepSize,
                                        sig = control$sig,
                                        gam = control$gam,
                                        maxIterOut = control$maxIterOut,
                                        maxIterIn = control$maxIterIn,
                                        maxIterLine = control$maxIterLine,
                                        epsOut = control$epsOut,
                                        epsIn = control$epsIn,
                                        useMultipleConvergenceCriteria = control$useMultipleConvergenceCriteria,
                                        regM2LLChangeEps = control$regM2LLChangeEps,
                                        verbose = control$verbose)
    
    SEM <- result$SEM
    
    parameterEstimates[i, names(parameters)] <- aCV4SEM:::getParameters(result$SEM, raw = FALSE)[names(parameters)]
    fits$m2LL[i] <- result$m2LL
    fits$regM2LL[i] <- result$regM2LL
    fits$convergence[i] <- result$convergence
    
    initialHessian4Optimizer <- result$Hessian - individualPenaltyFunctionHessian(aCV4SEM:::getParameters(result$SEM, raw = FALSE)[names(parameters)], 
                                                                                  currentTuningParameters,
                                                                                  penaltyFunctionArguments)
    
    if(control$saveHessian) Hessians$Hessian[[it]] <- initialHessian4Optimizer # note: we are returning the Hessian only for the likelihood, not for the penalty
    
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEMWithCustomPenalty",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(parameters),
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}

#' regularizeSEMWithCustomPenaltyRsolnp
#' 
#' optimize a SEM with custom penalty functions using Rsolnp optimizer (see ?Rsolnp::solnp). This optimizer is the default in regsem, for example
#' 
#' @param lavaanModel model of class lavaan 
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param tuningParameters data.frame with tuning parameter values. Important: The function will iterate over the rows of these tuning parameters and pass them to your penalty function
#' @param penaltyFunctionArguments arguments passed to individualPenaltyFunction, individualPenaltyFunctionGradient, and individualPenaltyFunctionHessian
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param saveHessian should the Hessian of the optimizer be saved? This will take a lot of space!
#' @param control option to set parameters of the optimizer; see ?Rsolnp::solnp
#' @export
regularizeSEMWithCustomPenaltyRsolnp <- function(lavaanModel, 
                                                 individualPenaltyFunction,
                                                 tuningParameters,
                                                 penaltyFunctionArguments,
                                                 startingValues = "est",
                                                 saveHessian = FALSE,
                                                 control = list("trace" = 0)){
  
  inputArguments <- as.list(environment())
  
  if(!is.data.frame(tuningParameters)) stop("tuningParameters must be a data frame (e.g., data.frame(lambda = seq(0,1,.1))).")
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  sampleSize <- nrow(rawData)
  
  ### initialize model ####
  if(any(startingValues == "est")){
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "est")
  }else if(any(startingValues == "start")){
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start")
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(aCV4SEM::getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See aCV4SEM::getLavaanParameters(lavaanModel).")
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start", 
                                   fit = FALSE)
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(startingValues), value = startingValues, raw = FALSE)
    SEM <- try(aCV4SEM:::fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to regularizeSEMWithCustomPenaltyRsolnp.")
  }
  
  # get parameters
  parameters <- aCV4SEM:::getParameters(SEM, raw = TRUE)
  
  # define fitting function
  fitfun <- function(parameters, 
                     SEM, 
                     sampleSize,
                     individualPenaltyFunction, 
                     tuningParameters,
                     penaltyFunctionArguments){
    SEM <- try(aCV4SEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters, raw = TRUE), silent = TRUE)
    if(is(SEM, "try-error")){
      return(99999999999999999)
    }
    SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)){
      return(99999999999999999)
    }
    
    penalty <- sampleSize*individualPenaltyFunction(parameters, tuningParameters, penaltyFunctionArguments)
    if(is(penalty, "try-error") || !is.finite(penalty)){
      return(99999999999999999)
    }
    
    return(SEM$m2LL + penalty)
  }
  
  fits <- data.frame(
    "m2LL" = rep(NA,nrow(tuningParameters)),
    "regM2LL"= rep(NA,nrow(tuningParameters)),
    "convergence" = rep(NA,nrow(tuningParameters))
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,nrow = nrow(tuningParameters), ncol = length(parameters)))
  colnames(parameterEstimates) <- names(parameters)
  parameterEstimates <- cbind(
    tuningParameters,
    parameterEstimates
  )
  
  if(saveHessian){
    Hessians <- list(
      "lambda" = tuningGrid$lambda,
      "alpha" = tuningGrid$alpha,
      "Hessian" = lapply(1:nrow(tuningGrid), 
                         matrix, 
                         data= NA, 
                         nrow=nrow(initialHessian), 
                         ncol=ncol(initialHessian))
    )
  }else{
    Hessians <- list(NULL)
  }
  
  progressbar = txtProgressBar(min = 0, 
                               max = nrow(tuningParameters), 
                               initial = 0, 
                               style = 3)
  
  parametersInit <- parameters
  
  for(i in 1:nrow(tuningParameters)){
    
    setTxtProgressBar(progressbar,i)
    
    currentTuningParameters <- tuningParameters[i,,drop = FALSE]
    
    result <- Rsolnp::solnp(pars = parametersInit, fun = fitfun, SEM = SEM, 
                            sampleSize = sampleSize,
                            individualPenaltyFunction = individualPenaltyFunction, 
                            tuningParameters = currentTuningParameters,
                            penaltyFunctionArguments = penaltyFunctionArguments,
                            control = control)
    parametersInit <- result$pars
    
    SEM <- try(aCV4SEM:::setParameters(SEM = SEM, labels = names(parametersInit), values = parametersInit, raw = TRUE), silent = TRUE)
    SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
    
    parameterEstimates[i, names(parameters)] <- aCV4SEM:::getParameters(SEM, raw = FALSE)[names(parameters)]
    fits$m2LL[i] <- SEM$m2LL
    fits$regM2LL[i] <- result$value[length(result$value)]
    fits$convergence[i] <- result$convergence == 0
    
    if(saveHessian) Hessians$Hessian[[it]] <- result$hessian # note: we are returning the Hessian for the likelihood AND the penalty
    
    if(result$convergence != 0) warning(paste0("Optimizer did not converge for ", paste0(names(currentTuningParameters), " = ",currentTuningParameters) , collapse = ", "))
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEMWithCustomPenalty",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(parameters),
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}