#### This script provides functions to perform approximate influence for different types of models ####

#' aI4regularizedSEM
#' 
#' Provides an approximate influence for models of class regularizedSEM. These models can be fit with regularizedSEM() (see ?aCV4SEM::regularizedSEM)
#' in this package.
#' 
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of subset fold. We recommend leave-one-out influence functions; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?aCV4SEM::aCV4regularizedSEM for an example
#' @param recomputeHessian if set to false, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed.
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGLMNET for more details
#' @examples
#' library(aCV4SEM)
#' 
#' # Let's first set up a regularized model. The following steps are
#' # explained in detail in ?aCV4SEM::regularizeSEM
#' dataset <- simulateExampleData()
#' 
#' lavaanSyntax <- "
#' f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
#'      l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
#'      l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
#' f ~~ 1*f
#' "
#' 
#' lavaanModel <- lavaan::sem(lavaanSyntax,
#'                            data = dataset,
#'                            meanstructure = TRUE,
#'                            std.lv = TRUE)
#' 
#' # Optional: Plot the model
#' # semPlot::semPaths(lavaanModel, 
#' #                   what = "est",
#' #                   fade = FALSE)
#' 
#' regsem <- regularizeSEM(
#'   lavaanModel = lavaanModel,
#'   regularizedParameterLabels = paste0("l", 6:15),
#'   penalty = "lasso",
#'   nLambdas = 5)
#' plot(regsem)
#' 
#' ## The approximate influence function can be computed with:
#' aI <- aI4regularizedSEM(regularizedSEM = regsem,
#'                           # we compute the influence of each individual
#'                           k = nrow(dataset)
#'                           )
#' 
#' # let's plot the influence
#' plot(aI) # The points are the fits for the model
#' # when one individual is removed. The lines are the
#' # fits with all individuals.
#' 
#' # To get more information, install the plotly - package and use::
#' # plot(aI, interactive  = TRUE)
#' # This creates an interactive plot which can be explored with your mouse
#' @export
aI4regularizedSEM <- function(regularizedSEM, k, recomputeHessian = TRUE, control = controlGLMNET()){
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  
  data <- try(lavaan::lavInspect(regularizedSEM@inputArguments$lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  N <- nrow(data)
  
  aCVSEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = regularizedSEM@inputArguments$lavaanModel, 
                                    transformVariances = TRUE, 
                                    fit = FALSE)
  
  # create subsets 
  if(is.matrix(k)){
    subsets <- k
    if(!is.logical(subsets)) stop("k must be a matrix with booleans (TRUE, FALSE)")
    if(nrow(subsets) != N) stop(paste0("k must have as many rows as there are subjects in your data set (", N, ")."))
    k <- ncol(subsets)
  }else{
    subsets <- aCV4SEM:::createSubsets(N = N, k = k)
  }
  
  # extract elements for easier access
  fits <- regularizedSEM@fits
  parameters <- regularizedSEM@parameters
  parameterLabels <- regularizedSEM@parameterLabels
  
  adaptiveLassoWeights <- regularizedSEM@inputArguments$adaptiveLassoWeights
  
  tuningParameters <- data.frame(lambda = fits$lambda, alpha = fits$alpha)
  
  regularizedParameterLabels <- regularizedSEM@inputArguments$regularizedParameterLabels
  
  subsetElements<- expand.grid(removedSubset = 1:k, 
                               lambda = unique(tuningParameters$lambda), 
                               alpha = unique(tuningParameters$alpha)
  )
  subsetFits <- cbind(subsetElements, 
                      matrix(NA, nrow(subsetElements), 
                             ncol = 1, 
                             dimnames = list(NULL, "fit")))
  
  subsetParameters <- cbind(subsetElements, 
                            matrix(NA, nrow(subsetElements), 
                                   ncol = length(parameterLabels), 
                                   dimnames = list(NULL, parameterLabels)))
  
  
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
      if(!regularizedSEM@inputArguments$control$saveHessian) stop("Hessians were not saved in the regularizedSEM object. This is the default as saving the Hessians can take a lot of disk space. To save the Hessians, use the controlGLMNET argument (see ?controlGLMNET).")
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
                                                                  hessianOfDifferentiablePart = hessianOfDifferentiablePart,
                                                                  control = control)
    subsetFits$fit[subsetFits$lambda == lambda & subsetFits$alpha == alpha] <- aInfluence$fits$fit
    subsetParameters[subsetFits$lambda == lambda & subsetFits$alpha == alpha, parameterLabels] <- aInfluence$subsetParameters[,parameterLabels]
    
    setTxtProgressBar(progressbar,ro)
    
  }
  
  return(
    new("approximateInfluence",
        subsets=subsets,
        subsetParameters = subsetParameters,
        tuningParameters = tuningParameters,
        subsetFits = subsetFits,
        parameters = regularizedSEM@parameters,
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
      subGroupGradient <- subGroupGradient + aCV4SEM::ridgeGradient(parameters = parameters,
                                                                    tuningParameters = penaltyFunctionTuning,
                                                                    penaltyFunctionArguments = penaltyFunctionArguments)
      if(is.null(hessianOfDifferentiablePart)){
        subGroupHessian <- subGroupHessian + aCV4SEM::ridgeHessian(parameters = parameters,
                                                                   tuningParameters = penaltyFunctionTuning,
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
                                       maxIter = control$maxIterIn, 
                                       epsBreak = control$epsIn, 
                                       useMultipleConvergenceCriteria = TRUE)
    
    rownames(direction) = names(parameters)
    
    # return step direction
    stepdirections[s,names(parameters)] <- direction[names(parameters),]
    
    # compute out of sample fit
    parameters_s <- parameters + stepdirections[s,names(parameters)]
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    SEM <- try(aCV4SEM:::fit(SEM), silent = TRUE)
    subsetParameters[subsetParameters$removedSubset == s, names(parameters)] <- aCV4SEM:::getParameters(SEM = SEM, raw = FALSE)[names(parameters)]
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
