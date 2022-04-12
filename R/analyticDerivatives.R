#' getScores
#' 
#' returns the scores of a model of class Rcpp_SEMCpp. This is the internal model representation. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' 
#' @param SEM model of class Rcpp_SEMCpp
#' @param raw controls if the internal transformations of aCV4SEM should be used. aCV4SEM will use an exponential function for all variances to 
#' avoid negative variances. When set to TRUE, the scores will be given for the internal parameter representation. Set to FALSE to get the usual 
#' scores
#' @export
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
#' @export
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
#' @export
getHessian <- function (SEM, raw = FALSE, eps = 1e-7){
  # THE FOLLOWING CODE IS ADAPTED FROM LAVAAN. 
  # SEE lavaan:::lav_model_hessian FOR THE IMPLEMENTATION
  # BY Yves Rosseel

  SEM <- fit(SEM = SEM)
  parameters <- getParameters(SEM = SEM, raw = raw)
  nParameters <- length(parameters)
  Hessian <- matrix(data = 0, 
                    nrow = nParameters, 
                    ncol = nParameters, 
                    dimnames = list(names(parameters),
                                    names(parameters)))
  
  for (p in 1:nParameters) {
    
    stepLeft <- twoStepLeft <- stepRight <- twoStepRight <- parameters
    stepLeft[p] <- stepLeft[p] - eps
    twoStepLeft[p] <- twoStepLeft[p] - 2 * eps
    stepRight[p] <- stepRight[p] + eps
    twoStepRight[p] <- twoStepRight[p] + 2 * eps
    
    SEM <- setParameters(SEM = SEM, labels = names(stepLeft), values = stepLeft, raw = raw)
    SEM <- fit(SEM = SEM)
    gradientsStepLeft <- getGradients(SEM = SEM, raw = raw)
    
    SEM <- setParameters(SEM = SEM, labels = names(twoStepLeft), values = twoStepLeft, raw = raw)
    SEM <- fit(SEM = SEM)
    gradientsTwoStepLeft <- getGradients(SEM = SEM, raw = raw)
    
    SEM <- setParameters(SEM = SEM, labels = names(stepRight), values = stepRight, raw = raw)
    SEM <- fit(SEM = SEM)
    gradientsStepRight <- getGradients(SEM = SEM, raw = raw)
    
    SEM <- setParameters(SEM = SEM, labels = names(twoStepRight), values = twoStepRight, raw = raw)
    SEM <- fit(SEM = SEM)
    gradientsTwoStepRight <- getGradients(SEM = SEM, raw = raw)
    
    Hessian[,p] <- (gradientsTwoStepLeft - 
                      8 * gradientsStepLeft + 
                      8 * gradientsStepRight - 
                      gradientsTwoStepRight)/(12 * eps)
  }
  # make symmetric
  Hessian <- (Hessian + t(Hessian))/2
  
  return(Hessian)
  
}