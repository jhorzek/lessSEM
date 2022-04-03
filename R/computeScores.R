computeScores <- function(SEM, individualFitfunction = individualMinus2LogLikelihood, raw, stepSize = 1e-5, ...){
  
  parameters <- getParameters(SEM, raw = raw)
  scores <- matrix(NA, nrow = nrow(SEM$data$rawData), ncol = length(parameters))
  colnames(scores) <- names(parameters)
  rownames(scores) <- paste0("person_", 1:nrow(SEM$data$rawData))
  
  for(parameter in 1:length(parameters)){
    
    # step forward
    parameters[parameter] <- parameters[parameter] + stepSize
    SEM <- setParameters(SEM, names(parameters), parameters, raw = raw)
    SEM <- computeExpected(SEM)
    for(i in 1:nrow(SEM$data$rawData)){
      scores[i,parameter] <- individualFitfunction(parameters, SEM, SEM$data$rawData[i,], raw = raw, ...)
    }
    
    # step backward
    parameters[parameter] <- parameters[parameter] - 2*stepSize
    SEM <- setParameters(SEM, names(parameters), parameters, raw = raw)
    SEM <- computeExpected(SEM)
    
    for(i in 1:nrow(SEM$data$rawData)){
      scores[i,parameter] <-  (scores[i,parameter] - individualFitfunction(parameters, SEM, SEM$data$rawData[i,], raw = raw, ...))/(2*stepSize)
    }
    
    # reset
    parameters[parameter] <- parameters[parameter] + stepSize
    
  }
  
  return(scores)
}
