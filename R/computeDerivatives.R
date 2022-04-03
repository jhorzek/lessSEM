computeGradients <- function(SEM, fitfunction = minus2LogLikelihood, stepSize = 1e-5, raw = FALSE, ...){
  
  parameters <- getParameters(SEM, raw = raw)
  gradients <- parameters
  gradients[] <- NA
  
  for(parameter in 1:length(parameters)){
    
    # step forward
    parameters[parameter] <- parameters[parameter] + stepSize
    SEM <- setParameters(SEM, names(parameters), parameters, raw = raw)
    SEM <- fit(SEM)
    gradients[parameter] <- fitfunction(parameters, SEM, raw = raw, ...)
    
    # step backward
    parameters[parameter] <- parameters[parameter] - 2*stepSize
    SEM <- setParameters(SEM, names(parameters), parameters, raw = raw)
    SEM <- fit(SEM)
    gradients[parameter] <- (gradients[parameter] - fitfunction(parameters, SEM, raw = raw, ...))/(2*stepSize)
    
    # reset
    parameters[parameter] <- parameters[parameter] + stepSize
    
  }
  
  return(gradients)
}

computeHessian <- function(SEM, fitfunction = minus2LogLikelihood, raw = FALSE, ...){
  
  parameters <- getParameters(SEM, raw = raw)
  
  hessian <- numDeriv::hessian(func = fitfunction, x = parameters, SEM = SEM, raw = raw, ... = ...)
  
  return(hessian)
}
