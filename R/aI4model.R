#### This script provides functions to perform approximate influence for different types of models ####

#' GLMNETApproximateInfluenceRcpp_SEMCpp
#' 
#' internal function for approximate influence based on the internal model representation of linr. The GLMNET part refers to the fact
#' that this function uses the GLMNET optimizer for non-differentiable penalty functions. 
#' 
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' @param subsets list with subsets created with createSubsets()
#' @param raw controls if the internal transformations of linr should be used.
#' @param regularizedParameterLabels vector with labels of regularized parameters
#' @param lambda value of tuning parameter lambda
#' @param alpha value of tuning parameter alpha (for elastic net)
#' @param adaptiveLassoWeights vector with labeled adaptive lasso weights. Only required if penalty = "adaptiveLasso"
#' @param maxIter maximal number of iterations used in GLMENT
#' @param epsBreak breaking condition of GLMNET
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGLMNET for more details
GLMNETApproximateInfluenceRcpp_SEMCpp <- function(SEM, 
                                                  subsets, 
                                                  raw = FALSE, 
                                                  regularizedParameterLabels,
                                                  lambda,
                                                  alpha = NULL,
                                                  adaptiveLassoWeights = NULL,
                                                  hessianOfDifferentiablePart = NULL,
                                                  control){
  
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  parameters <- linr:::getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  k <- ncol(subsets)
  
  # compute derivatives of -2log-Likelihood without penalty
  
  scores <- linr:::getScores(SEM = SEM, raw = raw)
  if(is.null(hessianOfDifferentiablePart)){
    hessian <- linr:::getHessian(SEM = SEM, raw = raw)
  }else{
    hessian <- hessianOfDifferentiablePart
  }
  
  # Now do the GLMNET inner iteration for each sub-group
  stepdirections <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(stepdirections) <- names(parameters)
  
  subsetParameters <- data.frame(lambda = rep(lambda,k), alpha = rep(alpha,k), removedSubset = 1:k)
  subsetParameters <- cbind(subsetParameters, 
                            matrix(NA, nrow = k, ncol = length(parameters), dimnames = list(NULL, names(parameters)))
  )
  
  fits <- data.frame(lambda = rep(lambda,k), alpha = rep(alpha,k), removedSubset = 1:k, fit = NA)
  
  for(s in 1:k){
    
    subGroupGradient <- apply(scores[-c(which(subsets[,s])),],2,sum) # gradients of all inidividuals but the ones in the subgroup
    subGroupHessian <- ((N-sum(subsets[,s]))/N)*hessian # approximated hessian for all but the subgroup
    subGroupLambda <- alpha*(N-sum(subsets[,s]))*lambda
    
    # if elastic net is used, we approximate the ridge part as well and add this here:
    if(alpha != 1){
      # ridge or elastic net are used
      if(length(unique(adaptiveLassoWeights))!= 1) warning("Combining ridge and elastic net with adaptive lasso weights is unexplored territory.")
      # we multiply with the adaptive lasso weights.
      penaltyFunctionTuning <- list("lambda" = unique(adaptiveLassoWeights)*lambda*(1-alpha)*(N-length(subsets[[s]])))
      penaltyFunctionArguments <- list("regularizedParameterLabels" = regularizedParameterLabels)
      subGroupGradient <- subGroupGradient + linr::ridgeGradient(parameters = parameters,
                                                                    tuningParameters = penaltyFunctionTuning,
                                                                    penaltyFunctionArguments = penaltyFunctionArguments)
      if(is.null(hessianOfDifferentiablePart)){
        subGroupHessian <- subGroupHessian + linr::ridgeHessian(parameters = parameters,
                                                                   tuningParameters = penaltyFunctionTuning,
                                                                   penaltyFunctionArguments = penaltyFunctionArguments)
      }
    }
    
    startDirection <- Sys.time()
    direction <- linr:::innerGLMNET(parameters = parameters, 
                                       N = (N-length(subsets[[s]])),
                                       subGroupGradient = subGroupGradient, 
                                       subGroupHessian = subGroupHessian, 
                                       subGroupLambda = subGroupLambda, 
                                       regularized = names(parameters)%in%regularizedParameterLabels, 
                                       adaptiveLassoWeights = adaptiveLassoWeights, 
                                       maxIter = control$maxIterIn, 
                                       epsBreak = control$epsIn, 
                                       useMultipleConvergenceCriteria = TRUE)
    
    rownames(direction) = names(parameters)
    
    # return step direction
    stepdirections[s,names(parameters)] <- direction[names(parameters),]
    
    # compute out of sample fit
    parameters_s <- parameters + stepdirections[s,names(parameters)]
    SEM <- linr:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    SEM <- try(linr:::fit(SEM), silent = TRUE)
    subsetParameters[subsetParameters$removedSubset == s, names(parameters)] <- linr:::getParameters(SEM = SEM, raw = FALSE)[names(parameters)]
    fits$fit[fits$removedSubset == s] <- SEM$m2LL 
    
    for(i in which(subsets[,s])){
      isMissing <- is.na(dataSet[i,])
      fits$fit[fits$removedSubset == s] <- fits$fit[fits$removedSubset == s] - computeIndividualM2LL(nObservedVariables = sum(!isMissing),  
                                                                                                     rawData = dataSet[i,!isMissing], 
                                                                                                     impliedMeans = SEM$impliedMeans[!isMissing], 
                                                                                                     impliedCovariance = SEM$impliedCovariance[!isMissing,!isMissing])
    }
  }
  
  return(
    list("fits" = fits,
         "subsets" = subsets,
         "scores" = scores,
         "subsetParameters" = subsetParameters)
  )
  
}
