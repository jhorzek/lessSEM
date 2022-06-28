#' CV4regularizedSEMWithCustomPenalty
#' 
#' Exact cross-validation for models of class regularizedSEMWithCustomPenalty. These models can be fit with regularizedSEMWithCustomPenalty.() (see ?regularizedSEMWithCustomPenalty.)
#' in this package.
#' 
#' @param regularizedSEMWithCustomPenalty. model of class regularizedSEMWithCustomPenalty.
#' @param k the number of cross-validation folds. Alternatively, 
#' a matrix with pre-defined subsets can be passed to the function. See ?linr::aCV4regularizedSEM for an example
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @examples 
#' library(linr)
#' 
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
#' tuningParameters <- data.frame(lambda = seq(0,1,.1)) # we will test 11 lambdas here
#' 
#' penaltyFunctionArguments <- list(
#'   eps = 1e-10,
#'   regularizedParameterLabels = paste0("l", 6:15)
#' )
#' 
#' smoothLASSO <- function(
#'   # here are our three arguments:
#'   parameters,
#'   tuningParameters,
#'   penaltyFunctionArguments
#' ){
#'   lambda <- tuningParameters$lambda 
#'   eps <- penaltyFunctionArguments$eps
#'   regularizedParameterLabels <- penaltyFunctionArguments$regularizedParameterLabels
#'   
#'   regularizedParameters <- parameters[regularizedParameterLabels]
#'   
#'   penaltyLasso <- lambda*sum(sqrt(regularizedParameters^2 + eps))
#'   
#'   return(penaltyLasso)
#' }
#' 
#' regsemApprox <- regularizeSEMWithCustomPenalty(lavaanModel = lavaanModel,
#'                                                individualPenaltyFunction = smoothLASSO,
#'                                                tuningParameters = tuningParameters,
#'                                                penaltyFunctionArguments = penaltyFunctionArguments)
#' 
#' ## true cross-validation
#' CV <- CV4regularizedSEMWithCustomPenalty(regularizedSEMWithCustomPenalty = regsemApprox,
#'                                          k = 5)
#' plot(CV)
#' @export
CV4regularizedSEMWithCustomPenalty <- function(regularizedSEMWithCustomPenalty, k, returnSubsetParameters = FALSE){
  if(!is(regularizedSEMWithCustomPenalty, "regularizedSEMWithCustomPenalty")){
    stop("regularizedSEMWithCustomPenalty must be of class regularizedSEMWithCustomPenalty")
  }
  
  data <- try(lavaan::lavInspect(regularizedSEMWithCustomPenalty@inputArguments$lavaanModel, "data"))
  if(is(data, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  N <- nrow(data)
  
  aCVSEM <- linr:::SEMFromLavaan(lavaanModel = regularizedSEMWithCustomPenalty@inputArguments$lavaanModel, transformVariances = TRUE, fit = FALSE)
  
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
  fits <- regularizedSEMWithCustomPenalty@fits
  parameters <- regularizedSEMWithCustomPenalty@parameters
  control <- regularizedSEMWithCustomPenalty@inputArguments$control
  
  # custom penalty bits
  individualPenaltyFunction <- regularizedSEMWithCustomPenalty@inputArguments$individualPenaltyFunction
  individualPenaltyFunctionGradient <- regularizedSEMWithCustomPenalty@inputArguments$individualPenaltyFunctionGradient
  individualPenaltyFunctionHessian <- regularizedSEMWithCustomPenalty@inputArguments$individualPenaltyFunctionHessian
  tuningParameters <- regularizedSEMWithCustomPenalty@inputArguments$tuningParameters
  penaltyFunctionArguments <- regularizedSEMWithCustomPenalty@inputArguments$penaltyFunctionArguments
  
  if(is.null(individualPenaltyFunctionGradient)){
    message("Using numDeriv to approximate the gradient of the individualPenaltyFunction.")
    individualPenaltyFunctionGradient <- linr::genericGradientApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  if(is.null(individualPenaltyFunctionHessian)){
    message("Using numDeriv to approximate the Hessian of the individualPenaltyFunction.")
    individualPenaltyFunctionHessian <- linr::genericHessianApproximiation
    penaltyFunctionArguments <- c(penaltyFunctionArguments, individualPenaltyFunction = individualPenaltyFunction)
  }
  
  cvfits <- data.frame(
    tuningParameters,
    cvfit = NA
  )
  
  cvfitsDetails <- as.data.frame(matrix(0,nrow = nrow(cvfits), ncol = k))
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
  
  # for CV - fits:
  SEM <- linr:::SEMFromLavaan(lavaanModel = regularizedSEMWithCustomPenalty@inputArguments$lavaanModel,
                                 whichPars = "est", 
                                 transformVariances = TRUE,
                                 fit = FALSE, 
                                 addMeans = TRUE, 
                                 activeSet = NULL
  )
  SEM <- linr:::fit(SEM)
  
  
  for(s in 1:k){
    cat("\n[",s, "/",k,"]\n")
    control_s <- control
    control_s$activeSet <- !subsets[,s]
    regularizedSEM_s <- linr::regularizeSEMWithCustomPenalty(
      lavaanModel = regularizedSEMWithCustomPenalty@inputArguments$lavaanModel,
      individualPenaltyFunction = individualPenaltyFunction,
      individualPenaltyFunctionGradient = individualPenaltyFunctionGradient,
      individualPenaltyFunctionHessian = individualPenaltyFunctionHessian, 
      tuningParameters = tuningParameters, 
      penaltyFunctionArguments = penaltyFunctionArguments,
      control = control_s
    )
    
    if(returnSubsetParameters){
      subsetParameters[s,,] <- regularizedSEM_s@parameters
    }
    
    for(p in 1:nrow(regularizedSEM_s@parameters)){
      
      for(i in which(subsets[,s])){
        
        cvfitsDetails[p, paste0("testSet",s)] <- cvfitsDetails[p, paste0("testSet",s)] + 
          linr:::individualMinus2LogLikelihood(par = unlist(regularizedSEM_s@parameters[p,regularizedSEM_s@parameterLabels]), 
                                                  SEM = SEM, 
                                                  data = data[i,,drop = FALSE], 
                                                  raw = FALSE)
      }
      
    }
    
  }
  
  cvfits$cvfit <- apply(cvfitsDetails[,paste0("testSet",1:k)],1,sum)
  
  return(
    new("CV4regularizedSEMWithCustomPenalty",
        parameters=parameters,
        tuningParameters = tuningParameters,
        cvfits = cvfits,
        parameterLabels = regularizedSEMWithCustomPenalty@parameterLabels,
        cvfitsDetails=cvfitsDetails, 
        subsets = subsets,
        subsetParameters = subsetParameters)
  )
}
