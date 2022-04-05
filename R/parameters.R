getParameters <- function(SEM, raw = FALSE){
  parameterTable <- SEM$getParameters()
  
  if(raw){
    values  <- parameterTable$rawValue
  }else{
    values  <- parameterTable$value
  }
  
  names(values) <- parameterTable$label
  
  values <- values[unique(names(values))]
  return(values)
}

setParameters <- function(SEM, labels, values, raw){
  if(length(labels) != length(values)){stop("labels and values of different length!")}
  SEM$setParameters(labels, values, raw)
  return(SEM)
}