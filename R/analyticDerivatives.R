#' .getScores
#' 
#' returns the scores of a model of class Rcpp_SEMCpp. This is the internal model 
#' representation. Models of this class
#' can be generated with the lessSEM:::.SEMFromLavaan-function.
#' 
#' @param SEM model of class Rcpp_SEMCpp
#' @param raw controls if the internal transformations of lessSEM should be used. 
#' lessSEM will use an exponential function for all variances to 
#' avoid negative variances. When set to TRUE, the scores will be given for the 
#' internal parameter representation. Set to FALSE to get the usual 
#' scores
#' @returns matrix with derivatives of the -2log-Likelihood for each person and parameter 
#' (rows are persons, columns are parameters)
.getScores <- function(SEM, raw){
  scores <- SEM$getScores(raw)
  colnames(scores) <- SEM$getParameterLabels()
  rownames(scores) <- paste0("person_",1:nrow(scores))
  return(scores)
}

#' .getGradients
#' 
#' returns the gradients of a model of class Rcpp_SEMCpp. This is the internal model 
#' representation. Models of this class
#' can be generated with the lessSEM:::.SEMFromLavaan-function.
#' 
#' @param SEM model of class Rcpp_SEMCpp
#' @param raw controls if the internal transformations of lessSEM should be used. 
#' lessSEM will use an exponential function for all variances to 
#' avoid negative variances. When set to TRUE, the gradients will be given for the 
#' internal parameter representation. Set to FALSE to get the usual 
#' gradients
#' @returns vector with derivatives of the -2log-Likelihood with respect to each parameter
.getGradients <- function(SEM, raw){
  gradients <- as.vector(SEM$getGradients(raw))
  names(gradients) <- names(.getParameters(SEM = SEM, raw = TRUE, transformations = FALSE))
  return(gradients)
}


#### Hessian ####
#' .getHessian
#' 
#' returns the Hessian of a model of class Rcpp_SEMCpp. This is the internal 
#' model representation. Models of this class
#' can be generated with the lessSEM:::.SEMFromLavaan-function. The function is adapted 
#' from lavaan::lav_model_hessian.
#' 
#' @param SEM model of class Rcpp_SEMCpp
#' @param raw controls if the internal transformations of lessSEM should be used. 
#' lessSEM will use an exponential function for all variances to 
#' avoid negative variances. When set to TRUE, the gradients will be given for 
#' the internal parameter representation. Set to FALSE to get the usual 
#' gradients
#' @param eps eps controls the step size of the numerical approximation.
#' @returns matrix with second derivatives of the -2log-Likelihood with respect 
#' to each parameter
.getHessian <- function (SEM, raw = FALSE, eps = 1e-7){
  # THE FOLLOWING CODE IS ADAPTED FROM LAVAAN. 
  # SEE lavaan:::lav_model_hessian FOR THE IMPLEMENTATION
  # BY Yves Rosseel
  
  SEM <- .fit(SEM = SEM)
  parameters <- .getParameters(SEM = SEM, raw = raw)
  hessian <- SEM$getHessian(names(parameters), parameters, raw, eps)
  rownames(hessian) <- names(parameters)
  colnames(hessian) <- names(parameters)
  return(hessian)
  
}
