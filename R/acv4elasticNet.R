#### This script provides functions to perform approximate cross-validation for different types of models ####

#' acv4lasso
#' 
#' approximate cross-validation for lasso regularized SEM
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?linr::aCV4regularizedSEM for an example
#' @param recomputeHessian if set to FALSE, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed. We currently recommend setting recomputeHessian to TRUE
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGlmnet for more details
#' @export
acv4lasso <- function(regularizedSEM, 
                      k, 
                      recomputeHessian = TRUE, 
                      returnSubsetParameters = FALSE,
                      control = controlGlmnet()
){
  return(
    acv4elasticNet(regularizedSEM = regularizedSEM, 
                   k = k, 
                   recomputeHessian = recomputeHessian, 
                   returnSubsetParameters = returnSubsetParameters,
                   control = control
    )
  )
}

#' acv4adaptiveLasso
#' 
#' approximate cross-validation for adaptive lasso regularized SEM
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?linr::aCV4regularizedSEM for an example
#' @param recomputeHessian if set to FALSE, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed. We currently recommend setting recomputeHessian to TRUE
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGlmnet for more details
#' @export
acv4adaptiveLasso <- function(regularizedSEM, 
                      k, 
                      recomputeHessian = TRUE, 
                      returnSubsetParameters = FALSE,
                      control = controlGlmnet()
){
  warning("Cross-validation will not recompute the adaptive lasso weights for the sub-samples! Make sure that this is appropriate.")
  return(
    acv4elasticNet(regularizedSEM = regularizedSEM, 
                   k = k, 
                   recomputeHessian = recomputeHessian, 
                   returnSubsetParameters = returnSubsetParameters,
                   control = control
    )
  )
}

#' acv4ridge
#' 
#' approximate cross-validation for ridge regularized SEM
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?linr::aCV4regularizedSEM for an example
#' @param recomputeHessian if set to FALSE, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed. We currently recommend setting recomputeHessian to TRUE
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGlmnet for more details
#' @export
acv4ridge <- function(regularizedSEM, 
                      k, 
                      recomputeHessian = TRUE, 
                      returnSubsetParameters = FALSE,
                      control = controlGlmnet()
){
  return(
    acv4elasticNet(regularizedSEM = regularizedSEM, 
                   k = k, 
                   recomputeHessian = recomputeHessian, 
                   returnSubsetParameters = returnSubsetParameters,
                   control = control
    )
  )
}

#' acv4elasticNet
#' 
#' approximate cross-validation for elastic net regularized SEM. These models can be fit with elasticNet() (see ?elasticNet)
#' in this package.
#' 
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?linr::aCV4regularizedSEM for an example
#' @param recomputeHessian if set to FALSE, the Hessians from the quasi newton optimization with GLMNET will be used. Otherwise the Hessian will be recomputed. We currently recommend setting recomputeHessian to TRUE
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @param control parameters passed to the GLMNET optimizer. Note that only arguments of the inner iteration are used. See ?controlGlmnet for more details
#' @examples
#' library(linr)
#' 
#' # Let's first set up a regularized model. The following steps are
#' # explained in detail in ?linr::regularizeSEM
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
#' ## The approximate cross-validation can be computed with:
#' aCV <- aCV4regularizedSEM(regularizedSEM = regsem,
#'                           # we highly recommend that you use approximate
#'                           # leave one out cross-validation:
#'                           k = nrow(dataset)
#' )
#' # let's plot the parameter values and the corresponding leave-one-out fit:
#' plot(aCV)
#' 
#' # the best parameters can be extracted with
#' coef(aCV)
#' 
#' # we can also use the one standard deviation rule:
#' coef(aCV, rule = "1sd")
#' 
#' # To see which person ended up in which sample, use:
#' aCV@subsets
#' 
#' #### Advanced ####
#' # If you are interested in the parameter estimates of each sub-sample,
#' # you must re-run the computation with returnSubsetParameters set to TRUE.
#' # This is disabled by default because it may take a lot of disk space
#' aCV <- aCV4regularizedSEM(regularizedSEM = regsem,
#'                           k = nrow(dataset),
#'                           returnSubsetParameters = TRUE
#' )
#' 
#' # the parameter are returned in a 3D-array. The subsets are in the
#' # rows, the parameters in the columns and the tuning parameters in
#' # the third dimension. To access the elements, use:
#' # Access elements for:
#' dimnames(aCV@subsetParameters)[[3]][1]
#' aCV@subsetParameters[,,1] # first lambda value
#' # Access elements for:
#' dimnames(aCV@subsetParameters)[[3]][2]
#' aCV@subsetParameters[,,2] # second lambda value
#' 
#' ## Currently, aCV4regularizedSEM recomputes the Hessian for each configuration
#' # of the tuning parameters. To reuse the Hessian of the previous optimiztation, use:
#' regsem <- regularizeSEM(
#'   lavaanModel = lavaanModel,
#'   regularizedParameterLabels = paste0("l", 6:15),
#'   penalty = "lasso",
#'   nLambdas = 5,
#'   control = controlGlmnet(saveHessian = TRUE)) # save the Hessian of the optimizer
#' 
#' aCV <- aCV4regularizedSEM(regularizedSEM = regsem,
#'                           k = nrow(dataset),
#'                           recomputeHessian = FALSE # don't recompute Hessian
#'                           
#' )
#' 
#' ## Instead of letting linr create the subsets, you can
#' # also pass your own subsets. As an example:
#' subsets <- linr:::createSubsets(N = nrow(dataset), k  = 10)
#' 
#' aCV <- aCV4regularizedSEM(regularizedSEM = regsem,
#'                           k = subsets
#' )
#' aCV@cvfitsDetails
#' @export
acv4elasticNet <- function(regularizedSEM, 
                           k, 
                           recomputeHessian = TRUE, 
                           returnSubsetParameters = FALSE,
                           control = controlGlmnet()
){
  
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  
  data <- try(lavaan::lavInspect(regularizedSEM@inputArguments$lavaanModel, "data"))
  if(is(data, "try-error")) 
    stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  N <- nrow(data)
  
  aCVSEM <- linr:::SEMFromLavaan(lavaanModel = regularizedSEM@inputArguments$lavaanModel, 
                                 transformVariances = TRUE, 
                                 fit = FALSE)
  
  # create subsets 
  if(is.matrix(k)){
    subsets <- k
    if(!is.logical(subsets)) stop("k must be a matrix with booleans (TRUE, FALSE)")
    if(nrow(subsets) != N) stop(paste0("k must have as many rows as there are subjects in your data set (", N, ")."))
    k <- ncol(subsets)
  }else{
    subsets <- linr:::createSubsets(N = N, k = k)
  }
  
  # extract elements for easier access
  fits <- regularizedSEM@fits
  parameters <- regularizedSEM@parameters
  
  weights <- regularizedSEM@inputArguments$weights
  
  tuningParameters <- data.frame(lambda = fits$lambda, 
                                 alpha = fits$alpha)
  
  cvfits <- data.frame(
    tuningParameters,
    cvfit = NA
  )
  
  cvfitsDetails <- as.data.frame(matrix(NA,
                                        nrow = nrow(cvfits), 
                                        ncol = k))
  colnames(cvfitsDetails) <- paste0("testSet", 1:k)
  cvfitsDetails <- cbind(
    tuningParameters,
    cvfitsDetails
  )
  
  if(returnSubsetParameters){
    subsetParameters <- array(NA, 
                              dim = c(k, length(regularizedSEM@parameterLabels), nrow(tuningParameters)),
                              dimnames = list(paste0("trainSet", 1:k),
                                              regularizedSEM@parameterLabels,
                                              NULL))
    dimname3 <- c()
    for(ro in 1:nrow(tuningParameters)){
      dimname3 <- c(dimname3, 
                    paste0(paste0(colnames(tuningParameters[ro,,drop = FALSE]),
                                  "=", 
                                  tuningParameters[ro,]), 
                           collapse = "; ")
      )
    }
    dimnames(subsetParameters)[[3]] <- dimname3
  }else{
    subsetParameters <- array(NA,dim = 1)
  }
  
  if(control$verbose == 0) progressbar = txtProgressBar(min = 0,
                                                        max = nrow(tuningParameters), 
                                                        initial = 0, 
                                                        style = 3)
  
  for(ro in 1:nrow(tuningParameters)){
    if(control$verbose != 0) cat("Model [",ro, "/", nrow(tuningParameters), "]\n")
    lambda <- tuningParameters$lambda[ro]
    alpha <- tuningParameters$alpha[ro]
    pars <- coef(regularizedSEM, lambda = lambda, alpha = alpha)
    
    aCVSEM <- linr:::setParameters(SEM = aCVSEM,
                                   labels = names(pars),
                                   values = pars,
                                   raw = FALSE)
    aCVSEM <- linr:::fit(aCVSEM)
    
    if(recomputeHessian){
      hessianOfDifferentiablePart <- NULL
    }else{
      if(!regularizedSEM@inputArguments$control$saveHessian)
        stop("Hessians were not saved in the regularizedSEM object. This is the default as saving the Hessians can take a lot of disk space. To save the Hessians, use the controlGlmnet argument (see ?controlGlmnet).")
      select <- regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$lambda == lambda &
        regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$alpha == alpha
      hessianOfDifferentiablePart <- regularizedSEM@internalOptimization$HessiansOfDifferentiablePart$Hessian[[which(select)]]
    }
    
    aCV <- try(linr:::acv4enet_GLMNET_SEMCpp(SEM = aCVSEM, 
                                             subsets = subsets,
                                             raw = TRUE, 
                                             weights = weights, 
                                             lambda = lambda,
                                             alpha = alpha,
                                             hessianOfDifferentiablePart = hessianOfDifferentiablePart,
                                             control = control))
    if(is(aCV, "try-error")) next
    
    if(returnSubsetParameters){
      subsetParameters[,,ro] <- aCV$subsetParameters[,dimnames(subsetParameters)[[2]]]
    }
    
    cvfitsDetails[cvfitsDetails$alpha == alpha & cvfitsDetails$lambda == lambda,
                  paste0("testSet",1:k)] <- aCV$leaveOutFits
    cvfits$cvfit[ro] <- sum(aCV$leaveOutFits)
    if(control$verbose == 0) setTxtProgressBar(progressbar,ro)
    
  }
  
  return(
    new("aCV4RegularizedSEM",
        parameters=parameters,
        cvfits = cvfits,
        parameterLabels = names(pars),
        regularized = names(weights)[weights != 0],
        cvfitsDetails=cvfitsDetails, 
        subsets = subsets,
        subsetParameters = subsetParameters)
  )
}

#' acv4enet_GLMNET_SEMCpp
#' 
#' internal function for approximate cross-validation based on the internal model representation of linr. The GLMNET part refers to the fact
#' that this function uses the GLMNET optimizer for non-differentiable penalty functions. 
#' 
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' @param subsets list with subsets created with createSubsets()
#' @param raw controls if the internal transformations of linr should be used.
#' @param weights labeled vector with weights of regularized parameters
#' @param lambda value of tuning parameter lambda
#' @param alpha value of tuning parameter alpha (for elastic net)
#' @param control arguments passed to GLMNET inner iteration
acv4enet_GLMNET_SEMCpp <- function(SEM, 
                                   subsets, 
                                   raw, 
                                   weights,
                                   lambda,
                                   alpha,
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
  if(control$verbose != 0) cat("Computing scores\n")
  scores <- linr:::getScores(SEM = SEM, raw = raw)
  if(is.null(hessianOfDifferentiablePart)){
    if(control$verbose != 0) cat("Computing Hessian\n\n")
    hessian <- linr:::getHessian(SEM = SEM, raw = raw)
    
    hessianEV <- eigen(hessian, only.values = TRUE)$values
    if(any(hessianEV < 0)) {
      hessian <- hessian + diag(-1.1*min(hessianEV), nrow(hessian))
      hessianEV <- eigen(hessian, only.values = TRUE)$values
      warning("Computed Hessian is not positive definite. Adding values to the diagonal to make it positive definite...")
      if(any(hessianEV < 0)) stop("Hessian still not positive definite...")
    }
    
  }else{
    hessian <- hessianOfDifferentiablePart
  }
  
  # Now do the GLMNET inner iteration for each sub-group
  stepdirections <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(stepdirections) <- names(parameters)
  
  subsetParameters <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(subsetParameters) <- names(parameters)
  rownames(subsetParameters) <- paste0("trainSet", 1:k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("testSet", 1:k)
  
  for(s in 1:k){
    Ntraining <- sum(!subsets[,s])
    Ntesting <- sum(subsets[,s])
    
    if(control$verbose != 0) cat("\r","Subset [",s, "/", k, "]")
    subGroupGradient <- apply(scores[-c(which(subsets[,s])),,drop = FALSE],2,sum) # gradients of all individuals but the ones in the subgroup
    subGroupHessian <- (Ntraining/N)*hessian # approximated hessian for all but the subgroup
    subGroupLambda <- alpha*(Ntraining)*lambda
    
    # if elastic net is used, we approximate the ridge part as well and add this here:
    if(alpha != 1){
      # ridge or elastic net are used
      if(length(unique(adaptiveLassoWeights))!= 1) warning("Combining ridge and elastic net with adaptive lasso weights is unexplored territory.")
      # we multiply with the adaptive lasso weights.
      penaltyFunctionTuning <- list("lambda" = weights[weights != 0]*lambda*(1-alpha)*(Ntraining))
      penaltyFunctionArguments <- list("regularizedParameterLabels" = names(weights)[weights != 0])
      subGroupGradient <- subGroupGradient + linr::ridgeGradient(parameters = parameters,
                                                                 tuningParameters = penaltyFunctionTuning,
                                                                 penaltyFunctionArguments = penaltyFunctionArguments)
      if(is.null(hessianOfDifferentiablePart)){
        subGroupHessian <- subGroupHessian + linr::ridgeHessian(parameters = parameters,
                                                                tuningParameters = penaltyFunctionTuning,
                                                                penaltyFunctionArguments = penaltyFunctionArguments)
      }
    }
    
    direction <- try(
      linr:::glmnetInner_C(parameters_kMinus1 = parameters, 
                           gradients_kMinus1 = subGroupGradient, 
                           Hessian = subGroupHessian, 
                           lambda = subGroupLambda, 
                           alpha = alpha, 
                           weights = weights, 
                           maxIterIn = Ntraining*control$maxIterIn,
                           breakInner = Ntraining*control$breakInner, 
                           verbose = control$verbose)
    )
    
    if(is(direction, "try-error")) next
    
    colnames(direction) = names(parameters)
    
    # return step direction
    stepdirections[s,names(parameters)] <- direction[,names(parameters)]
    
    # compute out of sample fit
    parameters_s <- parameters + stepdirections[s,names(parameters)]
    SEM <- linr:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    #SEM <- try(linr:::fit(SEM), silent = TRUE)
    subsetParameters[s,names(parameters)] <- linr:::getParameters(SEM = SEM, raw = FALSE)[names(parameters)]
    
    
    for(i in which(subsets[,s])){
      leaveOutFits[s] <- leaveOutFits[s] + linr:::individualMinus2LogLikelihood(par = subsetParameters[s,], 
                                                                                SEM = SEM, 
                                                                                data = dataSet[i,], 
                                                                                raw = FALSE)
    }
    
  }
  if(control$verbose != 0) cat("\n")
  return(
    list("leaveOutFits" = leaveOutFits,
         "subsets" = subsets,
         "scores" = scores,
         "subsetParameters" = subsetParameters)
  )
  
}