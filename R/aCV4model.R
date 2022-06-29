#### This script provides functions to perform approximate cross-validation for different types of models ####

#' aCV4lavaan
#' 
#' approximate cross-validation for models of class lavaan. 
#' 
#' @param lavaanModel model of class lavaan
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?lessSEM::aCV4regularizedSEM for an example
#' @param raw controls if the cross-validation should use the internal transformations of lessSEM. lessSEM will use an exponential function for all variances to 
#' avoid negative variances. This can result in better sub-group parameters, but may not be necessary and will also result in more difficult to interpret parameters.
#' @examples 
#' ## Approximate leave one out cross-validation for a lavaan model
#' ### set up model in lavaan
#' library(lavaan)
#' library(lessSEM)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'                     textual =~ x4 + x5 + x6
#'                     speed   =~ x7 + x8 + x9 '
#' HS <- HolzingerSwineford1939[,paste0("x",1:9)]
#' 
#' fit <- cfa(HS.model, data = HS, meanstructure = TRUE)
#' 
#' ### approximate cross-validation
#' aLOOCV <- aCV4lavaan(lavaanModel = fit, k = nrow(HS))
#' 
#' ### Optional: Compare this to a true leave one out cross-validation
#' # exactLOOCV <- rep(NA, nrow(HS))
#' # for(i in 1:nrow(HS)){
#' #   fit = sem(HS.model, HS[-i,], meanstructure = TRUE)
#' #   exactLOOCV[i] <- lessSEM:::computeIndividualM2LL(nObservedVariables = ncol(HS),
#' #                                                    rawData = as.numeric(HS[i,]),
#' #                                                    impliedMeans = fit@implied$mean[[1]],
#' #                                                    impliedCovariance = fit@implied$cov[[1]])
#' # }
#' # 
#' # # The plot shows the relation between exact and approximate cross-validation.
#' # # If the points are on the line, the approximate and exact cross-validation
#' # # produce relatively similar results
#' # plot(exactLOOCV, exactLOOCV, type = "l",
#' #      xlab = "exact loocv", ylab = "approximated loocv")
#' # points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
#' @export
aCV4lavaan <- function(lavaanModel, 
                       k, 
                       raw = FALSE){
  
  if(!is(lavaanModel, "lavaan")){
    stop("lavaanModel must be of class lavaan")
  }
  
  if(lavaanModel@Options$estimator != "ML") stop("lavaanModel must be fit with ml estimator.")
  
  N <- lavaan::lavInspect(lavaanModel, "nobs")
  aCVSEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, transformVariances = TRUE)
  
  # create subsets 
  if(is.matrix(k)){
    subsets <- k
    if(!is.logical(subsets)) stop("k must be a matrix with booleans (TRUE, FALSE)")
    if(nrow(subsets) != N) stop(paste0("k must have as many rows as there are subjects in your data set (", N, ")."))
    k <- ncol(subsets)
  }else{
    subsets <- lessSEM:::createSubsets(N = N, k = k)
  }
  
  individualPenaltyFunction <- function(a,b,c) {
    return(0)
  }
  individualPenaltyFunctionGradient <- function(a,b,c) {
    grs <- rep(0, length(a))
    names(grs) <- names(a)
    return(grs)
  }
  individualPenaltyFunctionHessian <- function(a,b,c) {
    return(
      matrix(0, 
             nrow = length(a), 
             ncol = length(a), 
             dimnames = list(names(a), names(a))))
  }
  
  aCV <- lessSEM:::customACVRcpp_SEMCpp(SEM = aCVSEM, 
                                        subsets = subsets,
                                        raw = raw, 
                                        individualPenaltyFunction = individualPenaltyFunction,
                                        individualPenaltyFunctionGradient = individualPenaltyFunctionGradient,
                                        individualPenaltyFunctionHessian = individualPenaltyFunctionHessian,
                                        currentTuningParameters = NULL,
                                        penaltyFunctionArguments = NULL, 
                                        hessianOfDifferentiablePart = NULL)
  return(aCV)
}

#' GLMNETACVRcpp_SEMCpp
#' 
#' internal function for approximate cross-validation based on the internal model representation of lessSEM. The GLMNET part refers to the fact
#' that this function uses the GLMNET optimizer for non-differentiable penalty functions. 
#' 
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' @param subsets list with subsets created with createSubsets()
#' @param raw controls if the internal transformations of lessSEM should be used.
#' @param regularizedParameterLabels vector with labels of regularized parameters
#' @param lambda value of tuning parameter lambda
#' @param alpha value of tuning parameter alpha (for elastic net)
#' @param adaptiveLassoWeights vector with labeled adaptive lasso weights. Only required if penalty = "adaptiveLasso"
#' @param control arguments passed to GLMNET inner iteration
GLMNETACVRcpp_SEMCpp <- function(SEM, 
                                 subsets, 
                                 raw, 
                                 regularizedParameterLabels,
                                 lambda,
                                 alpha = NULL,
                                 adaptiveLassoWeights = NULL,
                                 hessianOfDifferentiablePart = NULL,
                                 control){
  
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  parameters <- lessSEM:::getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  k <- ncol(subsets)
  
  # compute derivatives of -2log-Likelihood without penalty
  if(control$verbose != 0) cat("Computing scores\n")
  scores <- lessSEM:::getScores(SEM = SEM, raw = raw)
  if(is.null(hessianOfDifferentiablePart)){
    if(control$verbose != 0) cat("Computing Hessian\n\n")
    hessian <- lessSEM:::getHessian(SEM = SEM, raw = raw)
    
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
      penaltyFunctionTuning <- list("lambda" = unique(adaptiveLassoWeights)*lambda*(1-alpha)*(Ntraining))
      penaltyFunctionArguments <- list("regularizedParameterLabels" = regularizedParameterLabels)
      subGroupGradient <- subGroupGradient + lessSEM::ridgeGradient(parameters = parameters,
                                                                    tuningParameters = penaltyFunctionTuning,
                                                                    penaltyFunctionArguments = penaltyFunctionArguments)
      if(is.null(hessianOfDifferentiablePart)){
        subGroupHessian <- subGroupHessian + lessSEM::ridgeHessian(parameters = parameters,
                                                                   tuningParameters = penaltyFunctionTuning,
                                                                   penaltyFunctionArguments = penaltyFunctionArguments)
      }
    }
    
    direction <- try(lessSEM:::innerGLMNET(parameters = parameters, 
                                           N = Ntraining,
                                           subGroupGradient = subGroupGradient, 
                                           subGroupHessian = subGroupHessian, 
                                           subGroupLambda = subGroupLambda, 
                                           regularized = names(parameters)%in%regularizedParameterLabels, 
                                           adaptiveLassoWeights = adaptiveLassoWeights, 
                                           maxIter = control$maxIterIn, 
                                           epsBreak = control$epsIn, 
                                           useMultipleConvergenceCriteria = TRUE))
    if(is(direction, "try-error")) next
    
    rownames(direction) = names(parameters)
    
    # return step direction
    stepdirections[s,names(parameters)] <- direction[names(parameters),]
    
    # compute out of sample fit
    parameters_s <- parameters + stepdirections[s,names(parameters)]
    SEM <- lessSEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    #SEM <- try(lessSEM:::fit(SEM), silent = TRUE)
    subsetParameters[s,names(parameters)] <- lessSEM:::getParameters(SEM = SEM, raw = FALSE)[names(parameters)]
    
    
    for(i in which(subsets[,s])){
      leaveOutFits[s] <- leaveOutFits[s] + lessSEM:::individualMinus2LogLikelihood(par = subsetParameters[s,], 
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


#' aCV4regularizedSEMWithCustomPenalty
#' 
#' approximate cross-validation for models of class regularizedSEMWithCustomPenalty. These models can be fit with regularizedSEMWithCustomPenalty.() (see ?regularizedSEMWithCustomPenalty.)
#' in this package.
#' 
#' @param regularizedSEMWithCustomPenalty. model of class regularizedSEMWithCustomPenalty.
#' @param k the number of cross-validation folds. We recommend leave-one-out cross-validation; i.e. set k to the number of persons in the data set. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?lessSEM::aCV4regularizedSEM for an example
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @param recomputeHessian if set to FALSE, the Hessians from the quasi newton optimization with BFGS will be used. Otherwise the Hessian will be recomputed. We currently recommend setting recomputeHessian to TRUE
#' @export
aCV4regularizedSEMWithCustomPenalty <- function(regularizedSEMWithCustomPenalty,
                                                k, 
                                                returnSubsetParameters = FALSE,
                                                recomputeHessian = TRUE){
  if(!is(regularizedSEMWithCustomPenalty, "regularizedSEMWithCustomPenalty")){
    stop("regularizedSEMWithCustomPenalty must be of class regularizedSEMWithCustomPenalty")
  }
  
  data <- try(lavaan::lavInspect(regularizedSEMWithCustomPenalty@inputArguments$lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  N <- nrow(data)
  
  aCVSEM <- lessSEM:::SEMFromLavaan(lavaanModel = regularizedSEMWithCustomPenalty@inputArguments$lavaanModel, transformVariances = TRUE, fit = FALSE)
  
  # create subsets 
  if(is.matrix(k)){
    subsets <- k
    if(!is.logical(subsets)) stop("k must be a matrix with booleans (TRUE, FALSE)")
    if(nrow(subsets) != N) stop(paste0("k must have as many rows as there are subjects in your data set (", N, ")."))
    k <- ncol(subsets)
  }else{
    subsets <- lessSEM:::createSubsets(N = N, k = k)
  }
  
  # extract elements for easier access
  fits <- regularizedSEMWithCustomPenalty@fits
  parameters <- regularizedSEMWithCustomPenalty@parameters
  
  # custom penalty bits
  individualPenaltyFunction <- regularizedSEMWithCustomPenalty@inputArguments$individualPenaltyFunction
  individualPenaltyFunctionGradient <- regularizedSEMWithCustomPenalty@inputArguments$individualPenaltyFunctionGradient
  individualPenaltyFunctionHessian <- regularizedSEMWithCustomPenalty@inputArguments$individualPenaltyFunctionHessian
  tuningParameters <- regularizedSEMWithCustomPenalty@inputArguments$tuningParameters
  penaltyFunctionArguments <- regularizedSEMWithCustomPenalty@inputArguments$penaltyFunctionArguments
  
  if(is.null(individualPenaltyFunctionGradient)){
    message("Using numDeriv to approximate the gradient of the individualPenaltyFunction.")
    individualPenaltyFunctionGradient <- lessSEM::genericGradientApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  if(is.null(individualPenaltyFunctionHessian)){
    message("Using numDeriv to approximate the Hessian of the individualPenaltyFunction.")
    individualPenaltyFunctionHessian <- lessSEM::genericHessianApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  
  cvfits <- data.frame(
    tuningParameters,
    cvfit = NA
  )
  
  cvfitsDetails <- as.data.frame(matrix(NA,nrow = nrow(cvfits), ncol = k))
  colnames(cvfitsDetails) <- paste0("testSet", 1:k)
  cvfitsDetails <- cbind(
    tuningParameters,
    cvfitsDetails
  )
  
  if(returnSubsetParameters){
    subsetParameters <- array(NA, 
                              dim = c(k, length(regularizedSEMWithCustomPenalty@parameterLabels), nrow(tuningParameters)),
                              dimnames = list(paste0("trainSet", 1:k),
                                              regularizedSEMWithCustomPenalty@parameterLabels,
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
  
  progressbar = txtProgressBar(min = 0,
                               max = nrow(tuningParameters), 
                               initial = 0, 
                               style = 3)
  
  for(ro in 1:nrow(tuningParameters)){
    
    if(recomputeHessian){
      hessianOfDifferentiablePart <- NULL
    }else{
      
      if(is.null(regularizedSEMWithCustomPenalty@inputArguments$control$saveHessian) || !regularizedSEMWithCustomPenalty@inputArguments$control$saveHessian) stop("Hessians were not saved in the regularizedSEMWithCustomPenalty object. This is the default as saving the Hessians can take a lot of disk space. To save the Hessians, use the controlQuasiNewtonBFGS argument (see ?controlQuasiNewtonBFGS).")
      hessianOfDifferentiablePart <- regularizedSEMWithCustomPenalty@internalOptimization$HessiansOfDifferentiablePart$Hessian[[ro]]
      
    }
    
    currentTuningParameters <- tuningParameters[ro,,drop = FALSE]
    
    pars <- unlist(regularizedSEMWithCustomPenalty@parameters[ro,regularizedSEMWithCustomPenalty@parameterLabels])
    
    aCVSEM <- lessSEM:::setParameters(SEM = aCVSEM,
                                      labels = names(pars),
                                      values = pars,
                                      raw = FALSE)
    aCVSEM <- lessSEM:::fit(aCVSEM)
    
    aCV <- lessSEM:::customACVRcpp_SEMCpp(SEM = aCVSEM, 
                                          subsets = subsets,
                                          raw = TRUE, 
                                          individualPenaltyFunction = individualPenaltyFunction,
                                          individualPenaltyFunctionGradient = individualPenaltyFunctionGradient,
                                          individualPenaltyFunctionHessian = individualPenaltyFunctionHessian,
                                          currentTuningParameters = currentTuningParameters,
                                          penaltyFunctionArguments = penaltyFunctionArguments,
                                          hessianOfDifferentiablePart = hessianOfDifferentiablePart)
    
    if(returnSubsetParameters){
      subsetParameters[,,ro] <- aCV$subsetParameters[,dimnames(subsetParameters)[[2]]]
    }
    
    cvfitsDetails[ro,paste0("testSet",1:k)] <- aCV$leaveOutFits
    cvfits$cvfit[ro] <- sum(aCV$leaveOutFits)
    setTxtProgressBar(progressbar,ro)
    
  }
  
  return(
    new("aCV4regularizedSEMWithCustomPenalty",
        parameters=parameters,
        tuningParameters = tuningParameters,
        cvfits = cvfits,
        parameterLabels = names(pars),
        cvfitsDetails=cvfitsDetails, 
        subsets = subsets,
        subsetParameters = subsetParameters)
  )
}


#' customACVRcpp_SEMCpp
#' 
#' internal function for approximate cross-validation based on the internal model representation of lessSEM. The custom part refers to the fact
#' that this function uses a custom penalty function. 
#' 
#' @param SEM model of class Rcpp_SEMCpp. Models of this class
#' can be generated with the SEMFromLavaan-function.
#' @param subsets list with subsets created with createSubsets()
#' @param raw controls if the internal transformations of lessSEM should be used.
#' @param regularizedParameterLabels vector with labels of regularized parameters
#' @param lambda value of tuning parameter lambda
#' @param alpha value of tuning parameter alpha (for elastic net)
#' @param adaptiveLassoWeights vector with labeled adaptive lasso weights. Only required if penalty = "adaptiveLasso"
#' @param hessianOfDifferentiablePart hessian of log-likelihood
customACVRcpp_SEMCpp <- function(SEM, 
                                 subsets, 
                                 raw, 
                                 individualPenaltyFunction,
                                 individualPenaltyFunctionGradient,
                                 individualPenaltyFunctionHessian,
                                 currentTuningParameters,
                                 penaltyFunctionArguments,
                                 hessianOfDifferentiablePart){
  
  if(!is(SEM, "Rcpp_SEMCpp")){
    stop("SEM must be of class Rcpp_SEMCpp")
  }
  
  parameters <- lessSEM:::getParameters(SEM = SEM, raw = raw)
  dataSet <- SEM$rawData
  N <- nrow(dataSet)
  k <- ncol(subsets)
  
  # compute derivatives of -2log-Likelihood without penalty
  scores <- lessSEM:::getScores(SEM = SEM, raw = raw)
  if(is.null(hessianOfDifferentiablePart)){
    hessian <- lessSEM:::getHessian(SEM = SEM, raw = raw)
    
  }else{
    hessian <- hessianOfDifferentiablePart
  }
  
  # penalty scores
  penaltyScores <- individualPenaltyFunctionGradient(parameters, 
                                                     currentTuningParameters, 
                                                     penaltyFunctionArguments)
  penaltyHessian <- N*individualPenaltyFunctionHessian(parameters, 
                                                       currentTuningParameters, 
                                                       penaltyFunctionArguments)
  
  scores <- scores + matrix(rep(penaltyScores, N), 
                            nrow = N, 
                            ncol = length(penaltyScores), 
                            byrow = TRUE,
                            dimnames = list(NULL, names(penaltyScores)))[,colnames(scores)]
  hessian <- hessian + penaltyHessian[rownames(hessian), colnames(hessian)]
  
  hessianEV <- eigen(hessian, only.values = TRUE)$values
  if(any(hessianEV < 0)) {
    # see see https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
    hessian <- hessian + diag(-1.1*min(hessianEV), nrow(hessian))
    hessianEV <- eigen(hessian, only.values = TRUE)$values
    warning("Computed Hessian is not positive definite. Adding values to the diagonal to make it positive definite...")
    if(any(hessianEV < 0)) stop("Hessian still not positive definite...")
  }
  
  # Now do the optimization step for each sub-group
  stepdirections <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(stepdirections) <- names(parameters)
  
  subsetParameters <- matrix(NA, nrow = k, ncol = length(parameters))
  colnames(subsetParameters) <- names(parameters)
  rownames(subsetParameters) <- paste0("trainSet", 1:k)
  leaveOutFits <- rep(0, k)
  names(leaveOutFits) <- paste0("testSet", 1:k)
  
  for(s in 1:k){
    
    subGroupGradient <- apply(scores[-c(which(subsets[,s])),],2,sum) # gradients of all individuals but the ones in the subgroup
    subGroupHessian <- ((N-sum(subsets[,s]))/N)*hessian # approximated hessian for all but the subgroup
    
    direction <- -solve(subGroupHessian)%*%subGroupGradient
    
    rownames(direction) = names(parameters)
    
    # return step direction
    stepdirections[s,names(parameters)] <- direction[names(parameters),]
    
    # compute out of sample fit
    parameters_s <- parameters + stepdirections[s,names(parameters)]
    SEM <- lessSEM:::setParameters(SEM = SEM, labels = names(parameters), values = parameters_s, raw = raw)
    subsetParameters[s,names(parameters)] <- lessSEM:::getParameters(SEM = SEM, raw = FALSE)[names(parameters)]
    
    for(i in which(subsets[,s])){
      leaveOutFits[s] <- leaveOutFits[s] + lessSEM:::individualMinus2LogLikelihood(par = subsetParameters[s,], 
                                                                                   SEM = SEM, 
                                                                                   data = dataSet[i,], 
                                                                                   raw = FALSE)
    }
    
  }
  
  return(
    list("leaveOutFits" = leaveOutFits,
         "subsets" = subsets,
         "scores" = scores,
         "subsetParameters" = subsetParameters)
  )
  
}

#' createSubsets
#' 
#' create subsets for cross-validation
#' @param N number of samples in the data set
#' @param k number of subsets to create
#' @return matrix with subsets
createSubsets <- function(N,k){
  # build subgroups
  if(k < N){
    
    randomCases <- sample(1:N,N)
    subsets <- split(randomCases, sort(randomCases%%k)) # https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
    
  }else if(k == N){
    
    subsets <- vector("list",N)
    names(subsets) <- 1:N
    for(i in 1:N) subsets[[i]] <- i
    
  }else{
    stop(paste0("k must be <= ", N))
  }
  
  subsetMatrix <- matrix(NA, nrow = N, ncol = k,
                         dimnames = list(paste0("person",1:N), 
                                         paste0("testSet", 1:k)))
  for(s in 1:length(subsets)){
    subsetMatrix[,s] <- 1:N %in% subsets[[s]]
  }
  
  if(any(apply(subsetMatrix,1,sum) != 1)) stop("Error while splitting data in subsets: Some persons are in multiple or none of the subsets")
  
  return(subsetMatrix)
}