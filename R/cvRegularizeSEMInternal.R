#' .cvRegularizeSEMInternal
#' 
#' Combination of regularized structural equation model and cross-validation
#' 
#' Internal function: This function computes the regularized models
#' for all penalty functions which are implemented for glmnet and gist.
#' Use the dedicated penalty functions (e.g., lessSEM::cvLasso) to penalize
#' the model.
#' 
#' @param lavaanModel model of class lavaan 
#' @param k the number of cross-validation folds. Alternatively, a matrix with pre-defined subsets can be passed to the function. 
#' See ?lessSEM::cvLasso for an example
#' @param standardize should training and test sets be standardized?
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @param penalty string: name of the penalty used in the model
#' @param weights labeled vector with weights for each of the parameters in the 
#' model.
#' @param tuningParameters data.frame with tuning parameter values
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta() and controlGlmnet() functions.
#' @returns model of class cvRegularizedSEM
#' @import lavaan
#' @keywords internal
.cvRegularizeSEMInternal <- function(lavaanModel,
                                     k,
                                     standardize,
                                     penalty,
                                     weights,
                                     returnSubsetParameters,
                                     tuningParameters,
                                     method, 
                                     modifyModel,
                                     control){
  
  inputArguments <- as.list(environment())
  
  if(! method %in% c("ista", "glmnet")) 
    stop("Currently ony methods = 'ista' and methods = 'glmnet' are supported")
  if(method == "glmnet" & !penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet")) 
    stop(paste0(
      "glmnet only supports the following penalty functions: ",
      paste0(c("ridge", "lasso", "adaptiveLasso", "elasticNet"), collapse = ", ")
    )
    )
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet")
  
  if(method == "glmnet" && any(control$initialHessian == "lavaan")){
    rlang::inform(c("Note","Switching initialHessian from 'lavaan' to 'compute'."))
    control$initialHessian <- "compute"
  }
  
  if(!is(lavaanModel, "lavaan"))
    stop("lavaanModel must be of class lavaan")
  
  if(lavaanModel@Options$estimator != "ML") 
    stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) 
    stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  N <- nrow(rawData)
  
  misc <- list()
  
  # we have to set up the model once to find the parameters
  tmpSEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                           whichPars = "est",
                           addMeans = modifyModel$addMeans, 
                           activeSet = modifyModel$activeSet,
                           dataSet = modifyModel$dataSet,
                           transformations = modifyModel$transformations,
                           transformationList = modifyModel$transformationList,
                           transformationGradientStepSize = modifyModel$transformationGradientStepSize)
  
  
  parameterLabels <- names(.getParameters(SEM = tmpSEM, 
                                          raw = TRUE, 
                                          transformations = FALSE))
  
  # check if the data was standardized:
  if(all(apply(rawData, 2, function(x) abs(mean(x)) <= 1e-5))) 
    warning(paste0("It seems that you standardized your data before fitting your lavaan model. ",
                   "Note that this can undermine the results of the cross-validation due to dependencies between ",
                   "the subsets. Consider re-fitting your model with unstandardized data and use standardize = TRUE ",
                   "to automatically standardize the subsets.")
    )
  
  # create subsets 
  if(is.matrix(k)){
    subsets <- k
    if(!is.logical(subsets)) stop("k must be a matrix with booleans (TRUE, FALSE)")
    if(nrow(subsets) != N) stop(paste0("k must have as many rows as there are subjects in your data set (", N, ")."))
    k <- ncol(subsets)
  }else{
    subsets <- createSubsets(N = N, k = k)
  }
  
  if(penalty == "adaptiveLasso") 
    rlang::inform(c("Note",paste0("Automatic cross-validation for adaptiveLasso requested. ", 
                      "Note that using weights which are based on the full sample ",
                      "may undermine cross-validation. If the default is used (weights = NULL), ",
                      "weights for each subset will be computed using the inverse of the absolute MLE. ",
                      "Alternatively, pass a matrix as weights argument with weights for each subset.")
    ))
  
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
    # Holger Brandl at 
    # https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
    subsetParameters <- merge(tuningParameters,
                              data.frame(
                                trainSet = 1:k
                              ), 
                              by=NULL)
    
    subsetParameters <- cbind(
      subsetParameters,
      matrix(NA,
             nrow = nrow(subsetParameters),
             ncol = length(parameterLabels),
             dimnames = list(NULL, parameterLabels)
      )
    )
    
  }else{
    subsetParameters <- data.frame(NA)
  }
  
  if(penalty == "adaptiveLasso"){
    # save weights for inspection
    misc$newWeights <- data.frame(
      trainSet = 1:k)
    misc$newWeights <- cbind(misc$newWeights,
                             matrix(NA,
                                    nrow = k, 
                                    ncol = length(parameterLabels),
                                    dimnames = list(NULL, parameterLabels)
                             )
    )
  }
  
  for(s in 1:k){
    cat("\n[",s, "/",k,"]\n")
    control_s <- control
    # we need to pass the subset as our data set;
    # if scaling is used, this must also be applied here
    trainSet <- rawData[!subsets[,s],,drop = FALSE]
    testSet <- rawData[subsets[,s],,drop = FALSE]
    
    if(standardize){
      if(sum(subsets[,s]) < 2 || sum(!subsets[,s]) < 2){
        warning("Subsets too small for standardization Skipping standardization.")
      }else{
        # It is important to not scale the data prior to the splitting
        # Otherwise the data sets are not truly independent!
        rlang::inform(c("Note","Standardizing data sets ..."))
        trainSet <- scale(trainSet, center = TRUE, scale = TRUE)
        
        means <- attr(trainSet, "scaled:center")
        standardDeviations <- attr(trainSet, "scaled:scale")
        
        testSet <- cvScaler(testSet = testSet,
                            means = means, 
                            standardDeviations = standardDeviations)
      }
    }
    
    modifyModel$dataSet <- trainSet
    
    # check weights for adaptive Lasso
    ## option one: user provided weights
    if(penalty == "adaptiveLasso" && is.numeric(weights)){
      if(is.vector(weights)){
        warning("Using the same weights for all cv folds.")
        weights_s <- weights
      }
      if(is.matrix(weights) && nrow(weights) == k){
        weights_s <- weights[s,]
      }else{
        stop("Could not set the weights. Expected either a vector or a matrix with k rows")
      }
    }
    
    # option 2: the default
    if(penalty == "adaptiveLasso" && is.character(weights)){
      # weights are only specifying the names of the regularized parameters
      # use default;
      print(paste0("Computing new weights for sample ", s, "."))
      control_s$startingValues <- "start" # will result in new weights
      weights_s <- weights
      
    }
    if(penalty != "adaptiveLasso"){
      weights_s <- weights
    }
    
    regularizedSEM_s <- .regularizeSEMInternal(lavaanModel = lavaanModel, 
                                               penalty = penalty, 
                                               weights = weights_s, 
                                               tuningParameters = tuningParameters, 
                                               method = method,
                                               modifyModel = modifyModel,
                                               control = control_s
    )
    
    if(penalty == "adaptiveLasso"){
      
      misc$newWeights[misc$newWeights$trainSet == s,
                      names(regularizedSEM_s@weights)] <- regularizedSEM_s@weights
      
    }
    
    if(returnSubsetParameters){
      for(ro in 1:nrow(tuningParameters)){
        
        sel <- apply(subsetParameters[,colnames(tuningParameters)], 1, function(x) all(x == tuningParameters[ro,]))
        sel <- sel & subsetParameters$trainSet == s
        if(sum(sel) != 1) stop("Something went wrong while saving the subset parameters")
        subsetParameters[sel,parameterLabels] <- 
          as.matrix(regularizedSEM_s@parameters[ro,parameterLabels])
      }
    }
    
    # to compute the out of sample fit, we also need a SEM with all individuals 
    # if the test set
    SEM_s <- .SEMFromLavaan(
      lavaanModel = lavaanModel,
      whichPars = "start",
      fit = FALSE, 
      addMeans = modifyModel$addMeans,
      activeSet = modifyModel$activeSet, 
      dataSet = testSet,
      transformations = modifyModel$transformations,
      transformationList = modifyModel$transformationList
    )
    
    for(p in 1:nrow(regularizedSEM_s@parameters)){
      SEM_s <- .setParameters(
        SEM = SEM_s, 
        labels =  
          names(unlist(regularizedSEM_s@parameters[p,regularizedSEM_s@parameterLabels])),
        values = unlist(regularizedSEM_s@parameters[p,regularizedSEM_s@parameterLabels]),
        raw = FALSE)
      cvfitsDetails[p, paste0("testSet",s)] <- SEM_s$fit()
      
    }
    
  }
  
  cvfits$cvfit <- apply(cvfitsDetails[,paste0("testSet",1:k)],1,sum)
  
  # now, fit the full model
  
  tp <- tuningParameters[which.min(cvfits$cvfit)[1],]
  if(standardize) rawData <- scale(rawData)
  modifyModel$dataSet <- rawData
  regularizedSEM_full <- .regularizeSEMInternal(lavaanModel = lavaanModel, 
                                                penalty = penalty, 
                                                weights = weights_s, 
                                                tuningParameters = tp, 
                                                method = method,
                                                modifyModel = modifyModel,
                                                control = control
  )
  
  return(
    new("cvRegularizedSEM",
        parameters=regularizedSEM_full@parameters,
        transformations = regularizedSEM_full@transformations,
        cvfits = cvfits,
        parameterLabels = regularizedSEM_full@parameterLabels,
        regularized = regularizedSEM_full@regularized,
        cvfitsDetails = cvfitsDetails, 
        subsets = subsets,
        subsetParameters = subsetParameters,
        misc = misc)
  )
}


#' cvScaler
#' 
#' uses the means and standard deviations of the training set to standardize
#' the test set. See, e.g., https://scikit-learn.org/stable/modules/cross_validation.html .
#' 
#' @param testSet test data set
#' @param means means of the training set
#' @param standardDeviations standard deviations of the training set
#' @examples 
#' library(lessSEM)
#' data <- matrix(rnorm(50),10,5)
#' 
#' cvScaler(testSet = data, 
#'          means = 1:5, 
#'          standardDeviations = 1:5)
#' @returns scaled test set
#' @export
cvScaler <- function(testSet, means, standardDeviations){
  if(any(names(means) != colnames(testSet))) stop("Mismatch in names of means and testSet.")
  if(any(names(standardDeviations) != colnames(testSet))) stop("Mismatch in names of standardDeviations and testSet.")
  
  centered <- t(apply(testSet, 1, function(x) x-means))
  scaled <- t(apply(centered, 1, function(x) x/standardDeviations))
  colnames(scaled) <- colnames(testSet)
  
  return(scaled)
}