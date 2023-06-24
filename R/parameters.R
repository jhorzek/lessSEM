#' .getParameters
#' 
#' returns the parameters of the internal model representation.
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' @param raw controls if the parameter are returned in raw format or transformed
#' @param transformations should transformed parameters be included?
#' @returns labeled vector with parameter values
#' @keywords internal
.getParameters <- function(SEM, raw = FALSE, transformations = FALSE){
  
  if(is(SEM, "Rcpp_mgSEM")){
    parameterTable <- SEM$getParameters()
    if(transformations & raw){
      return(parameterTable$parmeters)
    }
    
    if(!transformations & raw){
      if(any(parameterTable$isTransformation)){
        return(parameterTable$parmeters[!parameterTable$isTransformation])
      }
      return(parameterTable$parmeters)
    }
    
    submodelParameters <- SEM$getSubmodelParameters()
    values <- c()
    
    if(transformations & !raw){
      for(sm in submodelParameters){
        addValues <- sm$value
        names(addValues) <- sm$label
        values  <- c(values, addValues)
      }
      values <- values[unique(names(values))] # remove duplicated elements
      
      # remove the raw elements from the parameter vector
      param <- c(parameterTable$parmeters,
                 values[!names(values) %in% names(parameterTable$parmeters)])
      # now replace the raw elements:
      param[names(values)] <- values
      
      return(param)
    }
    
    if(!transformations & !raw){
      for(sm in submodelParameters){
        addValues <- sm$value
        names(addValues) <- sm$label
        values  <- c(values, addValues)
      }
      values <- values[unique(names(values))] # remove duplicated elements
      
      # remove the raw elements from the parameter vector
      param <- c(parameterTable$parmeters,
                 values[!names(values) %in% names(parameterTable$parmeters)])
      # now replace the raw elements:
      param[names(values)] <- values
      
      # remove transformations
      if(any(parameterTable$isTransformation)){
        param <- param[names(param) %in%
                         names(parameterTable$parmeters)[!parameterTable$isTransformation]
        ]
      }
      
      return(param)
    }
    
  }
  
  if(is(SEM, "Rcpp_SEMCpp")){
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
  stop("Unknown model")
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