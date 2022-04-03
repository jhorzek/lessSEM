computeRegularizedDerivatives <- function(par, SEM, penaltyFunction, penaltyFunctionDerivative = NULL, ...){
  
  SEM <- setParameters(SEM = SEM, labels = names(par), values = par)
  parameters <- getParameters(SEM = SEM)
  
  # compute gradients of -2 log Likelihood
  gradientsM2LL <- computeAnalyticGradients(par = parameters, SEM = SEM)
  
  # compute gradients of penalty function
  if(is.null(penaltyFunctionDerivative)){
    
    penaltyGradients <- numericalDerivativeOfPenaltyFunction(parameters, penaltyFunction, ...)
    
  }else{
    
    penaltyGradients <- penaltyFunctionDerivative(parameters, ...)
    if(is.null(names(penaltyGradients))) stop("Can't match names of values returned by penaltyFunction and the gradients of the -2log-likelihood. Make sure that your penalty function returns the names of the parameters as well")
    
  }
  
  return(
    gradientsM2LL + penaltyGradients[names(gradientsM2LL)]
  )
  
}

computeRegularizedHessian <- function(par, SEM, raw, penaltyFunction, penaltyFunctionDerivative = NULL, ...){
  additionalArguments <- list(...)
  SEM <- setParameters(SEM = SEM, labels = names(par), values = par, raw = raw)
  parameters <- getParameters(SEM = SEM, raw = raw)
  
  HessianM2LL <- computeHessianFromAnalytic(par = parameters, SEM = SEM, raw = raw)
  
  # compute Hessian of penalty function
  if(!is.null(penaltyFunctionDerivative)){
    
    penaltyHessian <- do.call(what = numericalHessianOfPenaltyFunction, 
                              args = c(list("par" = parameters, 
                                            "raw" = raw, 
                                            "penaltyFunctionDerivative" = penaltyFunctionDerivative),
                                       additionalArguments
                              ))
    
  }else{
    
    penaltyHessian <- do.call(what = numericalHessianOfPenaltyFunction, 
                              args = c(list("par" = parameters, 
                                            "raw" = raw, 
                                            "penaltyFunctionDerivative" = aCV4SEM::numericalDerivativeOfPenaltyFunction,
                                            "penaltyFunction"  = penaltyFunction),
                                       additionalArguments
                              ))
    
    if(is.null(rownames(penaltyHessian))) stop("Can't match names of values returned by penaltyFunctionDerivative and the Hessian of the -2log-likelihood. Make sure that your penalty function returns the names of the parameters as well")
    
  }
  
  return(
    HessianM2LL + penaltyHessian[names(parameters), names(parameters)]
  )
  
}

## some helper functions
numericalDerivativeOfPenaltyFunction <- function(par, penaltyFunction, eps = 1e-7, ...){
  additionalArguments <- list(...)
  nParameters <- length(par)
  gradients <- rep(NA, nParameters)
  names(gradients) <- names(par)
  
  for(p in 1:nParameters){
    par_left <- par_right <- par
    par_left[p] <- par_left[p] - eps; par_right[p] <- par_right[p] + eps;
    gradients[p] <- (do.call(penaltyFunction, c(list("par" = par_right),
                                                additionalArguments)) - 
                       do.call(penaltyFunction, c(list("par" = par_left),
                                                  additionalArguments))) / (2*eps)
  }
  
  return(gradients)
}

numericalHessianOfPenaltyFunction <- function(par, raw, penaltyFunctionDerivative, eps = 1e-7, ...){
  additionalArguments <- list(...)
  
  nParameters <- length(par)
  Hessian <- matrix(data = 0, 
                    nrow = nParameters, 
                    ncol = nParameters, 
                    dimnames = list(names(par),
                                    names(par)))
  
  for (p in 1:nParameters) {
    
    stepLeft <- twoStepLeft <- stepRight <- twoStepRight <- par
    stepLeft[p] <- stepLeft[p] - eps
    twoStepLeft[p] <- twoStepLeft[p] - 2 * eps
    stepRight[p] <- stepRight[p] + eps
    twoStepRight[p] <- twoStepRight[p] + 2 * eps
    
    
    
    gradientsStepLeft <- do.call(penaltyFunctionDerivative, c(list("par" = stepLeft),
                                                              additionalArguments))
    gradientsTwoStepLeft <- do.call(penaltyFunctionDerivative, c(list("par" = twoStepLeft),
                                                                 additionalArguments))
    gradientsStepRight <- do.call(penaltyFunctionDerivative, c(list("par" = stepRight),
                                                               additionalArguments))
    gradientsTwoStepRight <- do.call(penaltyFunctionDerivative, c(list("par" = twoStepRight),
                                                                  additionalArguments))
    
    Hessian[,p] <- (gradientsTwoStepLeft - 
                      8 * gradientsStepLeft + 
                      8 * gradientsStepRight - 
                      gradientsTwoStepRight)/(12 * eps)
  }
  # make symmetric
  Hessian <- (Hessian + t(Hessian))/2
}