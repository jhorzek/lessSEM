#' getScores
#' 
#' returns the scores of a model of class Rcpp_SEMCpp. This is the internal model representation. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' 
#' @param SEM model of class Rcpp_SEMCpp
#' @param raw controls if the internal transformations of aCV4SEM should be used. aCV4SEM will use an exponential function for all variances to 
#' avoid negative variances. When set to TRUE, the scores will be given for the internal parameter representation. Set to FALSE to get the usual 
#' scores
getScores <- function(SEM, raw){
  scores <- SEM$getScores(raw)
  colnames(scores) <- SEM$getParameterLabels()
  rownames(scores) <- paste0("person_",1:nrow(scores))
  return(scores)
}

#' getGradients
#' 
#' returns the gradients of a model of class Rcpp_SEMCpp. This is the internal model representation. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' 
#' @param SEM model of class Rcpp_SEMCpp
#' @param raw controls if the internal transformations of aCV4SEM should be used. aCV4SEM will use an exponential function for all variances to 
#' avoid negative variances. When set to TRUE, the gradients will be given for the internal parameter representation. Set to FALSE to get the usual 
#' gradients
getGradients <- function(SEM, raw){
  gradients <- as.vector(SEM$getGradients(raw))
  names(gradients) <- SEM$getParameterLabels()
  return(gradients)
}


#### Hessian ####
#' getHessian
#' 
#' returns the Hessian of a model of class Rcpp_SEMCpp. This is the internal model representation. Models of this class
#' can be generated with the SEMFromLavaan-function. The function is adapted from \link[lavaan]{lav_model_hessian}.
#' 
#' @param SEM model of class Rcpp_SEMCpp
#' @param raw controls if the internal transformations of aCV4SEM should be used. aCV4SEM will use an exponential function for all variances to 
#' avoid negative variances. When set to TRUE, the gradients will be given for the internal parameter representation. Set to FALSE to get the usual 
#' gradients
#' @param eps eps controls the step size of the numerical approximation.
getHessian <- function (SEM, raw = FALSE, eps = 1e-7){
  # THE FOLLOWING CODE IS ADAPTED FROM LAVAAN. 
  # SEE lavaan:::lav_model_hessian FOR THE IMPLEMENTATION
  # BY Yves Rosseel
  
  SEM <- aCV4SEM:::fit(SEM = SEM)
  parameters <- aCV4SEM:::getParameters(SEM = SEM, raw = raw)
  hessian <- SEM$getHessian(names(parameters), parameters, raw, eps)
  rownames(hessian) <- names(parameters)
  colnames(hessian) <- names(parameters)
  return(hessian)
  
}
