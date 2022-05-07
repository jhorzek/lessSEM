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
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionHessian Hessian of the penalty function. Function should take the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' return a matrix with (length as parameters)^2 number of elements. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param tuningParameters data.frame with tuning parameter values. Important: The function will iterate over the rows of these tuning parameters and pass them to your penalty function
#' @param penaltyFunctionArguments arguments passed to individualPenaltyFunction, individualPenaltyFunctionGradient, and individualPenaltyFunctionHessian
#' @param method optimizer method. See ?optim
#' @param control control optimizer. See ?optim
#' @param hessian Should the optimizer return the hessian?
#' @param control option to set parameters of the optimizer
#' @export
regularizeSEMWithCustomPenalty <- function(lavaanModel, 
                                           individualPenaltyFunction,
                                           individualPenaltyFunctionGradient,
                                           individualPenaltyFunctionHessian,
                                           tuningParameters,
                                           penaltyFunctionArguments,
                                           control = aCV4SEM::controlQuasiNewtonBFGS()){
  
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
    "regM2LL"= rep(NA,nrow(tuningParameters))
  )
  
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