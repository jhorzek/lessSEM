#' .getParameters
#' 
#' returns the parameters of the internal model representation.
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' @param raw controls if the parameter are returned in raw format or transformed
#' @param transformations should transformed parameters be included?
#' @returns labeled vector with parameter values
#' @keywords internal
.getParameters <- function(SEM, raw = FALSE, transformations = FALSE){
  parameterTable <- SEM$getParameters()
  
  if(raw){
    values  <- parameterTable$rawValue
  }else{
    values  <- parameterTable$value
  }
  
  names(values) <- parameterTable$label
  
  values <- values[SEM$getParameterLabels()]
  
  if(!transformations) values <- values[names(values) %in% parameterTable$label[!parameterTable$isTransformation]]
  
  return(values)
}

#' .setParameters
#' 
#' change the parameters of the internal model representation.
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' @param labels vector with parameter labels
#' @param values vector with parameter values
#' @param raw are the parameters given in raw format or transformed?
#' @returns SEM with changed parameter values
#' @keywords internal
.setParameters <- function(SEM, labels, values, raw){
  if(length(labels) != length(values)){stop("labels and values of different length!")}
  SEM$setParameters(labels, values, raw)
  return(SEM)
}