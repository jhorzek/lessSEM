#' lsp
#' 
#' This function allows for regularization of models built in lavaan with the
#' lsp penalty.
#' The returned object is an S4 class; its elements can be accessed
#' with the "@" operator (see examples).
#' 
#' For more details, see:
#' 
#' 1. Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A general iterative shrinkage and thresholding algorithm for non-convex
#'  regularized optimization problems. Proceedings of the 30th International 
#'  Conference on Machine Learning, 28(2)(2), 37–45.
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta (see ?controlIsta)
#' @md
#' @examples 
#' library(lessSEM)
#' 
#' # Identical to regsem, lessSEM builds on the lavaan
#' # package for model specification. The first step
#' # therefore is to implement the model in lavaan.
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
#' # Optional: Plot the model
#' # semPlot::semPaths(lavaanModel,
#' #                   what = "est",
#' #                   fade = FALSE)
#' 
#' regsem <- lsp(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20),
#'   thetas = seq(0.01,2,length.out = 5))
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' 
#' # 5-fold cross-Validation
#' cv <- cv4regularizedSEM(regularizedSEM = regsem,
#'                           k = 5)
#' coef(cv)
#' @export
lsp <- function(lavaanModel,
                regularized,
                lambdas,
                thetas,
                modifyModel = lessSEM::modifyModel(),
                control = controlIsta()){
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- lessSEM::regularizeSEMInternal(lavaanModel = lavaanModel, 
                                 penalty = "lsp", 
                                 weights = regularized,
                                 tuningParameters = expand.grid(lambda = lambdas, 
                                                                theta = thetas), 
                                 method = "ista", 
                                 modifyModel = modifyModel, 
                                 control = control
  )
  
  return(result)
  
}

#' cv4lsp
#' 
#' Exact cross-validation for models of class regularizedSEM which have been fitted with lsp penalty
#' These models can be fit with lsp() (see ?lsp)
#' in this package.
#' 
#' @param regularizedSEM model of class regularizedSEM
#' @param k the number of cross-validation folds. Alternatively, a matrix with pre-defined subsets can be passed to the function. 
#' See ?lessSEM::aCV4regularizedSEM for an example
#' @param dataSet optional: Pass the full, unscaled data set to the function. 
#' This is important if the data has to be scaled prior to the analysis. If scaling
#' is performed on the full sample, this will result in dependencies between the 
#' subsets created by the cross-validation. To prevent this, pass the full data and use scaleData = TRUE,
#' the scalingFunction, and the scalingArguments
#' @param scaleData if set to TRUE, the subsets will be scaled using the scalingFunction
#' @param scalingFunction this function is used to scale the subsets. It MUST take two arguments:
#' first, the data set as matrix and second the scalingArguments. The latter can be anything you need
#' for the scaling
#' @param scalingArguments the second argument passed to scalingFunction. Can contain any number
#' of arguments needed for the scaling
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @examples
#' # see ?lsp
#' 
#' @export
cv4lsp <- function(regularizedSEM, 
                   k, 
                   dataSet = NULL,
                   scaleData = FALSE,
                   scalingFunction = function(dataSet,scalingArguments) 
                     scale(x = dataSet, 
                           center = scalingArguments$center, 
                           scale = scalingArguments$scale),
                   scalingArguments = list("center" = TRUE, 
                                           "scale" = TRUE),
                   returnSubsetParameters = FALSE){
  
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  
  lavaanData <- try(lavaan::lavInspect(regularizedSEM@inputArguments$lavaanModel, "data"))
  if(is(lavaanData, "try-error")) stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  
  if(is.null(dataSet)){
    message("Reusing the data set from regularizedSEM. If your data was rescaled prior to the analysis, this may result in incorrect cross-validation results. Pass the data using the dataset argument and the scalingFunction to let this function rescale your data for each subset.")
    data <- lavaanData
  }else{
    data <- dataSet[,colnames(lavaanData),drop = FALSE]
    if(any(!colnames(lavaanData) %in% colnames(data))) 
      stop("Not all variables present in the lavaanModel are in the dataSet")
  }
  
  N <- nrow(data)
  
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
  fits <- regularizedSEM@fits
  parameters <- regularizedSEM@parameters
  control <- regularizedSEM@inputArguments$control
  weights <- regularizedSEM@inputArguments$weights
  modifyModel <- regularizedSEM@inputArguments$modifyModel
  if(any(!weights %in% c(0,1))) stop("Weights must be 0 or 1.")
  
  tuningParameters <- data.frame(theta = fits$theta,
                                 lambda = fits$lambda
  )
  
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
  
  
  for(s in 1:k){
    cat("\n[",s, "/",k,"]\n")
    control_s <- control
    # we need to pass the subset as our data set;
    # if scaling is used, this must also be applied here
    trainSet <- data[!subsets[,s],,drop = FALSE]
    testSet <- data[subsets[,s],,drop = FALSE]
    
    if(scaleData){
      if(sum(subsets[,s]) < 2 || sum(!subsets[,s]) < 2){
        warning("Subsets too small for scaling. Skipping scaling.")
      }else{
        # It is important to not scale the data prior to the splitting
        # Otherwise the data sets are not truly independent!
        message("Scaling data sets ...")
        trainSet <- scalingFunction(trainSet, scalingArguments)
        testSet <- scalingFunction(testSet, scalingArguments)
      }
    }
    
    modifyModel$dataSet <- trainSet
    
    weights_s <- regularizedSEM@inputArguments$weights
    
    regularizedSEM_s <- lessSEM::lsp(
      lavaanModel = regularizedSEM@inputArguments$lavaanModel,
      regularized = names(weights_s)[weights_s!=0],
      lambdas = regularizedSEM@inputArguments$lambdas, 
      thetas = regularizedSEM@inputArguments$thetas,
      modifyModel = modifyModel,
      control = control_s
    )
    
    if(returnSubsetParameters){
      for(ro in 1:nrow(tuningParameters)){
        dimname3 <- paste0(paste0(colnames(tuningParameters[ro,,drop = FALSE]),
                                  "=", 
                                  tuningParameters[ro,]), 
                           collapse = "; ")
        
        subsetParameters[s,,dimname3] <- as.matrix(regularizedSEM_s@parameters[ro,dimnames(subsetParameters)[[2]]])
      }
    }
    
    # to compute the out of sample fit, we also need a SEM with all individuals 
    # if the test set
    SEM_s <- lessSEM::SEMFromLavaan(
      lavaanModel = regularizedSEM@inputArguments$lavaanModel,
      whichPars = "start",
      fit = FALSE, 
      addMeans = modifyModel$addMeans,
      activeSet = modifyModel$activeSet, 
      dataSet = testSet, 
      transformVariances = TRUE
    )
    
    for(p in 1:nrow(regularizedSEM_s@parameters)){
      SEM_s <- lessSEM::setParameters(
        SEM = SEM_s, 
        labels =  
          names(unlist(regularizedSEM_s@parameters[p,regularizedSEM_s@parameterLabels])),
        values = unlist(regularizedSEM_s@parameters[p,regularizedSEM_s@parameterLabels]),
        raw = FALSE)
      cvfitsDetails[p, paste0("testSet",s)] <- SEM_s$fit()
      
      # for(i in which(subsets[,s])){
      #   
      #   cvfitsDetails[p, paste0("testSet",s)] <- cvfitsDetails[p, paste0("testSet",s)] + 
      #     lessSEM:::individualMinus2LogLikelihood(par = unlist(regularizedSEM_s@parameters[p,regularizedSEM_s@parameterLabels]), 
      #                                             SEM = SEM, 
      #                                             data = data[i,,drop = FALSE], 
      #                                             raw = FALSE)
      # }
      
    }
    
  }
  
  cvfits$cvfit <- apply(cvfitsDetails[,paste0("testSet",1:k)],1,sum)
  
  return(
    new("CV4RegularizedSEM",
        parameters=parameters,
        cvfits = cvfits,
        parameterLabels = regularizedSEM@parameterLabels,
        regularized = regularizedSEM@parameterLabels[regularizedSEM@inputArguments$weights != 0],
        cvfitsDetails = cvfitsDetails, 
        subsets = subsets,
        subsetParameters = subsetParameters)
  )
}