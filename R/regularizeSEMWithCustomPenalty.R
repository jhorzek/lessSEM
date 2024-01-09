#' .regularizeSEMWithCustomPenaltyRsolnp
#' 
#' Optimize a SEM with custom penalty function using the Rsolnp optimizer (see ?Rsolnp::solnp). This optimizer is the default in regsem (see ?regsem::cv_regsem).
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
#' @param carryOverParameters should parameters from the previous iteration be used as starting values of
#' the next iteration?
#' @param control option to set parameters of the optimizer; see ?Rsolnp::solnp
#' @returns Model of class regularizedSEMWithCustomPenalty
#' @keywords internal
.regularizeSEMWithCustomPenaltyRsolnp <- function(lavaanModel, 
                                                  individualPenaltyFunction,
                                                  tuningParameters,
                                                  penaltyFunctionArguments,
                                                  startingValues = "est",
                                                  carryOverParameters = TRUE,
                                                  control = list("trace" = 0)){
  if(!("Rsolnp" %in% rownames(utils::installed.packages())))
    stop(".regularizeSEMWithCustomPenaltyRsolnp requires the Rsolnp package to be installed.")
  
  inputArguments <- as.list(environment())
  
  if(!is.data.frame(tuningParameters)) stop("tuningParameters must be a data frame (e.g., data.frame(lambda = seq(0,1,.1))).")
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(.getRawData(lavaanModel, NULL, "fiml")$rawData)
  if(is(rawData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  sampleSize <- nrow(rawData)
  
  ### initialize model ####
  if(any(startingValues == "est")){
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "est",
                          addMeans = control$addMeans, 
                          activeSet = control$activeSet)
  }else if(any(startingValues == "start")){
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "start",
                          addMeans = control$addMeans, 
                          activeSet = control$activeSet)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See lessSEM::getLavaanParameters(lavaanModel).")
    SEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                          whichPars = "start", 
                          fit = FALSE,
                          addMeans = control$addMeans, 
                          activeSet = control$activeSet)
    SEM <- .setParameters(SEM = SEM, 
                          labels = names(startingValues), 
                          values = startingValues, 
                          raw = FALSE)
    SEM <- try(.fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$objectiveValue)) stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to regularizeSEMWithCustomPenaltyRsolnp.")
  }
  
  # check if we have a likelihood objective function:
  usesLikelihood <- all(SEM$getEstimator() == "fiml")
  
  # get parameters
  parameters <- .getParameters(SEM, raw = TRUE)
  
  # define fitting function
  fitfun <- function(parameters, 
                     SEM, 
                     sampleSize,
                     individualPenaltyFunction, 
                     tuningParameters,
                     penaltyFunctionArguments){
    SEM <- try(.setParameters(SEM = SEM, 
                              labels = names(parameters), 
                              values = parameters, 
                              raw = TRUE), 
               silent = TRUE)
    if(is(SEM, "try-error")){
      return(99999999999999999)
    }
    SEM <- try(.fit(SEM), silent = TRUE)
    if(is(SEM, "try-error") || !is.finite(SEM$objectiveValue)){
      return(99999999999999999)
    }
    
    # check if covariance matrix is positive definite
    if(any(eigen(SEM$S, only.values = TRUE)$values < 0)){
      return(99999999999999999)
    }
    if(any(eigen(SEM$impliedCovariance, only.values = TRUE)$values < 0)){
      return(99999999999999999)
    }
    
    penalty <- sampleSize*individualPenaltyFunction(parameters, 
                                                    tuningParameters, 
                                                    penaltyFunctionArguments)
    if(is(penalty, "try-error") || !is.finite(penalty)){
      return(99999999999999999)
    }
    
    # division by N to be closer to the implementation in regsem
    return((SEM$objectiveValue + penalty)/sampleSize)
  }
  
  fits <- data.frame(
    "objectiveValue" = NA,
    "regObjectiveValue" = NA,
    "m2LL" = rep(NA,nrow(tuningParameters)),
    "regM2LL"= rep(NA,nrow(tuningParameters)),
    "convergence" = rep(NA,nrow(tuningParameters))
  )
  fits <- cbind(tuningParameters,
                fits)
  parameterEstimates <- as.data.frame(matrix(NA,nrow = nrow(tuningParameters), ncol = length(parameters)))
  colnames(parameterEstimates) <- names(parameters)
  parameterEstimates <- cbind(
    tuningParameters,
    parameterEstimates
  )
  
  Hessians <- list(NULL)
  
  progressbar = utils::txtProgressBar(min = 0, 
                                      max = nrow(tuningParameters), 
                                      initial = 0, 
                                      style = 3)
  
  parametersInit <- parameters
  
  for(i in 1:nrow(tuningParameters)){
    
    utils::setTxtProgressBar(progressbar,i)
    
    currentTuningParameters <- tuningParameters[i,,drop = FALSE]
    
    result <- try(Rsolnp::solnp(pars = parametersInit, fun = fitfun, SEM = SEM, 
                                sampleSize = sampleSize,
                                individualPenaltyFunction = individualPenaltyFunction, 
                                tuningParameters = currentTuningParameters,
                                penaltyFunctionArguments = penaltyFunctionArguments,
                                control = control))
    if(is(result, "try-error")) next
    # regsem does not re-use the parameters for the next iteration;
    # this is considerably slower, but it does help the optimizer.
    # The optimizer builds an approximation of the Hessian which
    # will be better if more iterations have to be done until convergence
    if(carryOverParameters) parametersInit <- result$pars
    
    SEM <- try(.setParameters(SEM = SEM, labels = names(result$pars), 
                              values = result$pars, raw = TRUE), 
               silent = TRUE)
    SEM <- try(.fit(SEM), silent = TRUE)
    
    parameterEstimates[i, names(parameters)] <- .getParameters(SEM, raw = FALSE)[names(parameters)]
    
    fits$objectiveValue[i] <- SEM$objectiveValue
    fits$regObjectiveValue[i] <- fits$objectiveValue[i] + sampleSize*individualPenaltyFunction(result$pars, 
                                                                                                 currentTuningParameters, 
                                                                                                 penaltyFunctionArguments)
    
    if(usesLikelihood){
      fits$m2LL[i] <- fits$objectiveValue[i]
      fits$regM2LL[i] <- fits$regObjectiveValue[i]
    }
    
    #result$value[length(result$value)]
    fits$convergence[i] <- result$convergence == 0
    
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
                 inputArguments = inputArguments,
                 notes = c("Notes:"))
  
  return(results)
  
}
