#### This script provides functions to perform approximate cross-validation for different types of models ####

#' aCV4regularizedSEM
#' 
#' approximate cross-validation for models of class regularizedSEM. These models can be fit with regularizedSEM() (see ?regularizedSEM)
#' in this package.
#' 
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param recomputeHessian if set to false, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed.
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGLMNET for more details
#' @export
aCV4regularizedSEM <- function(regularizedSEM, k, recomputeHessian = TRUE, control = controlGLMNET()){
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
  
  adaptiveLassoWeights <- regularizedSEM@inputArguments$adaptiveLassoWeights
  
  tuningParameters <- data.frame(lambda = fits$lambda, alpha = fits$alpha)
  
  regularizedParameterLabels <- regularizedSEM@inputArguments$regularizedParameterLabels
  
  cvfits <- data.frame(
    tuningParameters,
    cvfit = NA
  )
  
  cvfitsDetails <- as.data.frame(matrix(NA,nrow = nrow(cvfits), ncol = k))
  colnames(cvfitsDetails) <- paste0("sample", 1:k)
  cvfitsDetails <- cbind(
    tuningParameters,
    cvfitsDetails
  )
  
  subsets <- aCV4SEM:::createSubsets(N = N, k = k)
  
  if(control$verbose == 0) progressbar = txtProgressBar(min = 0,
                                                        max = nrow(tuningParameters), 
                                                        initial = 0, 
                                                        style = 3)
  
  for(ro in 1:nrow(tuningParameters)){
    if(control$verbose != 0) cat("Model [",ro, "/", nrow(tuningParameters), "]\n")
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
    
    aCV <- aCV4SEM:::GLMNETACVRcpp_SEMCpp(SEM = aCVSEM, 
                                          subsets = subsets,
                                          raw = TRUE, 
                                          regularizedParameterLabels = regularizedParameterLabels, 
                                          lambda = lambda,
                                          alpha = alpha, 
                                          adaptiveLassoWeights = adaptiveLassoWeights,
                                          hessianOfDifferentiablePart = hessianOfDifferentiablePart,
                                          control = control)
    
    cvfitsDetails[cvfitsDetails$alpha == alpha & cvfitsDetails$lambda == lambda,paste0("sample",1:k)] <- aCV$leaveOutFits
    cvfits$cvfit[ro] <- sum(aCV$leaveOutFits)
    if(control$verbose == 0) setTxtProgressBar(progressbar,ro)
    
  }
  
  return(
    new("aCV4RegularizedSEM",
        parameters=parameters,
        cvfits = cvfits,
        parameterLabels = names(pars),
        regularized = regularizedParameterLabels,
        cvfitsDetails=cvfitsDetails, 
        subsets = subsets)
  )
}

#' aCV4lavaan
#' 
#' approximate cross-validation for models of class lavaan. 
#' 
#' @param lavaanModel model of class lavaan
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set
#' @param raw controls if the cross-validation should use the internal transformations of aCV4SEM. aCV4SEM will use an exponential function for all variances to 
#' avoid negative variances. This can result in better sub-group parameters, but may not be necessary and will also result in more difficult to interpret parameters.
#' @export
aCV4lavaan <- function(lavaanModel, 
                       k, 
                       raw = FALSE){
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  aCVSEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = lavaanModel, transformVariances = TRUE)
  
  return(aCV4SEM:::smoothACVRcpp_SEMCpp(SEM = aCVSEM, 
                                        k = k,
                                        individualPenaltyFunction = NULL, 
                                        individualPenaltyFunctionGradient = NULL,
                                        individualPenaltyFunctionHessian = NULL,
                                        raw = raw, 
                                        penaltyFunctionArguments = NULL))
}