#### This script provides functions to perform approximate influence for different types of models ####

#' aI4regularizedSEM
#' 
#' approximate influence for models of class regularizedSEM. These models can be fit with regularizedSEM() (see ?regularizedSEM)
#' in this package.
#' 
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param recomputeHessian if set to false, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed.
#' @export
aI4regularizedSEM <- function(regularizedSEM, k, recomputeHessian = TRUE){
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  
  data <- try(lavaan::lavInspect(regularizedSEM@inputArguments$lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  N <- nrow(data)
  
  aCVSEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = regularizedSEM@inputArguments$lavaanModel, transformVariances = TRUE, fit = FALSE)
  
  # extract elements for easier access
  fits <- regularizedSEM@fits
  parameters <- regularizedSEM@parameters
  parameterLabels <- regularizedSEM@parameterLabels
  
  adaptiveLassoWeights <- regularizedSEM@inputArguments$adaptiveLassoWeights
  
  tuningParameters <- data.frame(lambda = fits$lambda, alpha = fits$alpha)
  
  regularizedParameterLabels <- regularizedSEM@inputArguments$regularizedParameterLabels
  
  subsets <- aCV4SEM:::createSubsets(N = N, k = k)
  
  subsetParameters <- array(NA, dim = c(k,length(parameterLabels), nrow(tuningParameters)), dimnames = list(paste0("subset",1:k), 
                                                                                                            parameterLabels, 
                                                                                                            NULL))
  subsetFits <- array(NA, dim = c(k,1, nrow(tuningParameters)), dimnames = list(paste0("subset",1:k), 
                                                                                "fit", 
                                                                                NULL))
  
  progressbar = txtProgressBar(min = 0,
                               max = nrow(tuningParameters), 
                               initial = 0, 
                               style = 3)
  
  for(ro in 1:nrow(tuningParameters)){
    
    lambda <- tuningParameters$lambda[ro]
    alpha <- tuningParameters$alpha[ro]
    pars <- coef(regularizedSEM, lambda = lambda, alpha = alpha)
    
    aCVSEM <- aCV4SEM:::setParameters(SEM = aCVSEM,
                                      labels = names(pars),
                                      values = pars,
                                      raw = FALSE)
    aCVSEM <- aCV4SEM:::fit(aCVSEM)
    
    if(recomputeHessian){
      hessianOfDifferentiablePart <- NULL
    }else{
      select <- regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$lambda == lambda &
        regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$alpha == alpha
      hessianOfDifferentiablePart <- regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$Hessian[[which(select)]]
    }
    
    aInfluence <- aCV4SEM:::GLMNETApproximateInfluenceRcpp_SEMCpp(SEM = aCVSEM, 
                                                                  subsets = subsets,
                                                                  raw = TRUE, 
                                                                  regularizedParameterLabels = regularizedParameterLabels, 
                                                                  lambda = lambda,
                                                                  alpha = alpha, 
                                                                  adaptiveLassoWeights = adaptiveLassoWeights,
                                                                  hessianOfDifferentiablePart = hessianOfDifferentiablePart)
    subsetFits[,,ro] <- aInfluence$fits
    subsetParameters[,,ro] <- aInfluence$subsetParameters[,parameterLabels]
    setTxtProgressBar(progressbar,ro)
    
  }
  
  return(
    new("approximateInfluence",
        subsets=subsets,
        subsetParameters = subsetParameters,
        tuningParameters = tuningParameters,
        subsetFits = subsetFits,
        parameterLabels = parameterLabels,
        regularized = regularizedParameterLabels)
  )
}


#' GLMNETApproximateInfluenceRcpp_SEMCpp
#' 
#' internal function for approximate influence based on the internal model representation of aCV4SEM. The GLMNET part refers to the fact
#' that this function uses the GLMNET optimizer for non-differentiable penalty functions. 
#' 
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' @param subsets list with subsets created with createSubsets()
#' @param raw controls if the internal transformations of aCV4SEM should be used.
#' @param regularizedParameterLabels vector with labels of regularized parameters
#' @param lambda value of tuning parameter lambda
#' @param alpha value of tuning parameter alpha (for elastic net)
#' @param adaptiveLassoWeights vector with labeled adaptive lasso weights. Only required if penalty = "adaptiveLasso"
#' @param maxIter maximal number of iterations used in GLMENT
#' @param epsBreak breaking condition of GLMNET
GLMNETApproximateInfluenceRcpp_SEMCpp <- function(SEM, 
                                                  subsets, 
                                                  raw = FALSE, 
                                                  regularizedParameterLabels,
                                                  lambda,
                                                  alpha = NULL,
                                                  adaptiveLassoWeights = NULL,
                                                  hessianOfDifferentiablePart = NULL,
                                                  maxIter = 100,
                                                  epsBreak = 1e-10){
  
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  parameters <- aCV4SEM:::getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  k <- ncol(subsets)
  
  # compute derivatives of -2log-Likelihood without penalty
  
  scores <- aCV4SEM:::getScores(SEM = SEM, raw = raw)
  if(is.null(hessianOfDifferentiablePart)){
    hessian <- aCV4SEM:::getHessian(SEM = SEM, raw = raw)
  }else{
    hessian <- hessianOfDifferentiablePart
  }
  
  # Now do the GLMNET inner iteration for each sub-group
  stepdirections <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(stepdirections) <- names(parameters)
  
  #subsetParameters <- vector("list",k)
  subsetParameters <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(subsetParameters) <- names(parameters)
  rownames(subsetParameters) <- paste0("sample", 1:k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("sample", 1:k)
  fits <- rep(NA, k)
  names(fits) <- paste0("FitWithoutSample", 1:k)
  
  for(s in 1:k){
    
    subGroupGradient <- apply(scores[-c(which(subsets[,s])),],2,sum) # gradients of all inidividuals but the ones in the subgroup
    subGroupHessian <- ((N-sum(subsets[,s]))/N)*hessian # approximated hessian for all but the subgroup
    subGroupLambda <- alpha*(N-sum(subsets[,s]))*lambda
    
    # if elastic net is used, we approximate the ridge part as well and add this here:
    if(alpha != 1){
      # ridge or elastic net are used
      if(length(unique(adaptiveLassoWeights))!= 1) warning("Combining ridge and elastic net with adaptive lasso weights is unexplored territory.")
      # we multiply with the adaptive lasso weights.
      penaltyFunctionArguments <- list("regularizedParameterLabels" = regularizedParameterLabels,
                                       "lambda" = unique(adaptiveLassoWeights)*lambda*(1-alpha)*(N-length(subsets[[s]])))
      subGroupGradient <- subGroupGradient + aCV4SEM::ridgeGradient(parameters = parameters,
                                                                    penaltyFunctionArguments = penaltyFunctionArguments)
      if(is.null(hessianOfDifferentiablePart)){
        subGroupHessian <- subGroupHessian + aCV4SEM::ridgeHessian(parameters = parameters,
                                                                   penaltyFunctionArguments = penaltyFunctionArguments)
      }
    }
    
    startDirection <- Sys.time()
    direction <- aCV4SEM:::innerGLMNET(parameters = parameters, 
                                       subGroupGradient = subGroupGradient, 
                                       subGroupHessian = subGroupHessian, 
                                       subGroupLambda = subGroupLambda, 
                                       regularized = names(parameters)%in%regularizedParameterLabels, 
                                       adaptiveLassoWeights = adaptiveLassoWeights, 
                                       maxIter = maxIter, 
                                       epsBreak = epsBreak)
    
    rownames(direction) = names(parameters)
    
    # return step direction
    stepdirections[s,names(parameters)] <- direction[names(parameters),]
    
    # compute out of sample fit
    parameters_s <- parameters + stepdirections[s,names(parameters)]
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
    subsetParameters[s,names(parameters)] <- aCV4SEM:::getParameters(SEM = SEM, raw = FALSE)[names(parameters)]
    fits[s] <- SEM$m2LL 
    
    for(i in which(subsets[,s])){
      isMissing <- is.na(dataSet[i,])
      fits[s] <- fits[s] - computeIndividualM2LL(nObservedVariables = sum(!isMissing),  
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
