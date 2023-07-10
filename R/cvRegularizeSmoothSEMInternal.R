#' .cvRegularizeSmoothSEMInternal
#' 
#' Combination of smoothly regularized structural equation model and cross-validation
#' 
#' Internal function: This function computes the regularized models
#' for all penalty functions which are implemented for bfgs.
#' Use the dedicated penalty functions (e.g., lessSEM::cvSmoothLasso) to penalize
#' the model.
#' 
#' @param lavaanModel model of class lavaan 
#' @param k the number of cross-validation folds. Alternatively, a matrix with pre-defined subsets can be passed to the function. 
#' See ?lessSEM::cvSmoothLasso for an example
#' @param standardize should training and test sets be standardized?
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @param penalty string: name of the penalty used in the model
#' @param weights labeled vector with weights for each of the parameters in the 
#' model.
#' @param tuningParameters data.frame with tuning parameter values
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method optimizer used. Currently only "bfgs" is supported.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @returns model of class cvRegularizedSEM
#' @keywords internal
.cvRegularizeSmoothSEMInternal <- function(lavaanModel,
                                           k,
                                           standardize,
                                           penalty,
                                           weights,
                                           returnSubsetParameters,
                                           tuningParameters,
                                           epsilon,
                                           modifyModel,
                                           method = "bfgs",
                                           control){
  
  notes <- c("Notes:")
  
  inputArguments <- as.list(environment())
  
  if(!penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet")) 
    stop(paste0(
      "bfgs only supports the following penalty functions: ",
      paste0(c("ridge", "lasso", "adaptiveLasso", "elasticNet"), collapse = ", ")
    )
    )
  
  if(!is(control, "controlBFGS")) 
    stop("control must be of class controlBfgs See ?controlBfgs.")
  
  if(!is(lavaanModel, "lavaan"))
    stop("lavaanModel must be of class lavaan")
  
  if(standardize & (!tolower(lavaanModel@Options$estimator) %in% c("ml", "fiml","mlm", "mlmv", "mlmvs", "mlf", "mlr")))
    stop("Automatic standardization is currently only implemented for maximum likelihood estimation.")
  
  control$breakOuter <- .adaptBreakingForWls(lavaanModel = lavaanModel, 
                                             currentBreaking = control$breakOuter,
                                             selectedDefault = ifelse(method == "ista",
                                                                      control$breakOuter == controlIsta()$breakOuter,
                                                                      control$breakOuter == controlGlmnet()$breakOuter
                                             ))
  
  rawData <- try(.getRawData(lavaanModel, NULL, "fiml")$rawData)
  if(is(rawData, "try-error")) 
    stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  N <- nrow(rawData)
  
  misc <- list()
  
  parameterLabels <- names(getLavaanParameters(lavaanModel))
  
  if(all(apply(rawData, 2, function(x) abs(mean(x, na.rm = TRUE)) <= 1e-5))) 
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
    notes <- c(
      notes,
      paste0("Automatic cross-validation for adaptiveLasso requested. ", 
             "Note that using weights which are based on the full sample ",
             "may undermine cross-validation. If the default is used (weights = NULL), ",
             "weights for each subset will be computed using the inverse of the absolute MLE. ",
             "Alternatively, pass a matrix as weights argument with weights for each subset.")
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
        notes <- c(
          notes, 
          "Using automatic standardizing for cross-validation data sets"
        )
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
    
    regularizedSEM_s <- .regularizeSmoothSEMInternal(lavaanModel = lavaanModel, 
                                                     penalty = penalty, 
                                                     weights = weights_s, 
                                                     tuningParameters = tuningParameters, 
                                                     epsilon = epsilon,
                                                     tau = 0,
                                                     modifyModel = modifyModel,
                                                     method = method,
                                                     control = control_s
    )
    notes <- c(notes,
               regularizedSEM_s@notes
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
      transformationList = modifyModel$transformationList,
      transformationGradientStepSize = modifyModel$transformationGradientStepSize
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
  regularizedSEM_full <- .regularizeSmoothSEMInternal(lavaanModel = lavaanModel, 
                                                      penalty = penalty, 
                                                      weights = weights_s, 
                                                      tuningParameters = tp, 
                                                      epsilon = epsilon,
                                                      tau = 0,
                                                      modifyModel = modifyModel,
                                                      control = control
  )
  
  notes <- unique(notes)
  
  if(length(notes) > 1){
    cat("\n")
    rlang::inform(notes)
  }
  
  return(
    new("cvRegularizedSEM",
        parameters=regularizedSEM_full@parameters,
        cvfits = cvfits,
        parameterLabels = regularizedSEM_full@parameterLabels,
        regularized = regularizedSEM_full@regularized,
        cvfitsDetails = cvfitsDetails, 
        subsets = subsets,
        subsetParameters = subsetParameters,
        misc = misc,
        notes = notes)
  )
}