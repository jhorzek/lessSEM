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
  
  notes <- c("Notes:")
  
  if(! method %in% c("ista", "glmnet")) 
    stop("Currently ony methods = 'ista' and methods = 'glmnet' are supported")
  if(method == "glmnet" & !penalty %in% c("ridge", "lasso", "adaptiveLasso", "elasticNet", "scad", "cappedL1", "mcp", "lsp")) 
    stop(paste0(
      "glmnet only supports the following penalty functions: ",
      paste0(c("ridge", "lasso", "adaptiveLasso", "elasticNet", "scad", "cappedL1", "mcp", "lsp"), collapse = ", ")
    )
    )
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet")
  
  if(method == "glmnet" && any(control$initialHessian == "lavaan")){
    notes <- c(notes, "Switching initialHessian from 'lavaan' to 'compute'.")
    control$initialHessian <- "compute"
  }
  
  if(!is(lavaanModel, "lavaan"))
    stop("lavaanModel must be of class lavaan")
  
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
  
  # we have to set up the model once to find the parameters
  tmpSEM <- .SEMFromLavaan(lavaanModel = lavaanModel,
                           whichPars = "est",
                           addMeans = modifyModel$addMeans, 
                           activeSet = modifyModel$activeSet,
                           dataSet = modifyModel$dataSet,
                           transformations = modifyModel$transformations,
                           transformationList = modifyModel$transformationList,
                           transformationGradientStepSize = modifyModel$transformationGradientStepSize)
  
  misc$estimator <- tmpSEM$getEstimator()
  parameterLabels <- names(.getParameters(SEM = tmpSEM, 
                                          raw = TRUE, 
                                          transformations = FALSE))
  
  # check if the data was standardized:
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
    notes <- c(notes, 
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
      
      if(!tolower(lavaanModel@Options$estimator) %in% c("ml", "fiml","mlm", "mlmv", "mlmvs", "mlf", "mlr"))
        stop("Automatic standardization is currently only implemented for maximum likelihood estimation.")
      
      if(sum(subsets[,s]) < 2 || sum(!subsets[,s]) < 2){
        warning("Subsets too small for standardization Skipping standardization.")
      }else{
        # It is important to not scale the data prior to the splitting
        # Otherwise the data sets are not truly independent!
        notes <- c(notes, 
                   "Using automatic standardization for cross-validation data sets."
                   )
        trainSet <- scale(trainSet, center = TRUE, scale = TRUE)
        
        means <- attr(trainSet, "scaled:center")
        standardDeviations <- attr(trainSet, "scaled:scale")
        
        testSet <- cvScaler(testSet = testSet,
                            means = means, 
                            standardDeviations = standardDeviations)
      }
    }
    
    # We update the lavaan models to set up our regularized model
    if(tolower(lavaanModel@Options$estimator) %in% c("ml", "fiml","mlm", "mlmv", "mlmvs", "mlf", "mlr")){
      
      # Note: lavaan will throw an error if N < p and missing != "ml", which we don't really care about
      # in case of ml-estimation because N does not have to be larger than p. In fact,
      # we can even use N = 1.
      # See e.g., p. 334 in 
      # Voelkle, M. C., Oud, J. H. L., von Oertzen, T., & Lindenberger, U. (2012).
      # Maximum Likelihood Dynamic Factor Modeling for Arbitrary N and T Using SEM. 
      # Structural Equation Modeling: A Multidisciplinary Journal, 19(3), 329–350. 
      # https://doi.org/10.1080/10705511.2012.687656
      
      lavaanModelTrain <- lavaanModel
      lavaanModelTest <- lavaanModel
      
      # replace data
      lavaanModelTrain@Data@X[[1]] <- as.matrix(trainSet)
      lavaanModelTrain@Options$do.fit <- FALSE # set to FALSE because otherwise lessSEM will
      # try to compare the fit of the model using only the train set to the fit of the
      # lavaanModelTrain, which was not fitted.
      
      lavaanModelTest@Data@X[[1]] <- as.matrix(testSet)
      lavaanModelTest@Options$do.fit <- FALSE
      
    }else{
      # If any other estimator is used, we need lavaan to set up the weight
      # matrices, etc for the WLS
      lavaanModelTrain <- lavaanModelTest <- lavaanModel
      lavaanModelTrain@Options$do.fit <- FALSE
      lavaanModelTest@Options$do.fit <- FALSE
      lavaanModelTrain <- suppressWarnings(.updateLavaan(lavaanModel = lavaanModel, 
                                                          key = "data",
                                                          value = trainSet))
      lavaanModelTest <- suppressWarnings(.updateLavaan(lavaanModel = lavaanModel, 
                                                         key = "data",
                                                         value = testSet))
    }
    
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
    
    regularizedSEM_s <- .regularizeSEMInternal(lavaanModel = lavaanModelTrain, 
                                               penalty = penalty, 
                                               weights = weights_s, 
                                               tuningParameters = tuningParameters, 
                                               method = method,
                                               modifyModel = modifyModel,
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
      lavaanModel = lavaanModelTest,
      whichPars = "start",
      fit = FALSE, 
      addMeans = modifyModel$addMeans,
      activeSet = modifyModel$activeSet, 
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
  if(tolower(lavaanModel@Options$estimator) %in% c("ml", "fiml","mlm", "mlmv", "mlmvs", "mlf", "mlr"))
  {
    if(standardize) rawData <- scale(rawData)
    modifyModel$dataSet <- rawData
  }else{
    modifyModel$dataSet <- NULL
  }
  regularizedSEM_full <- .regularizeSEMInternal(lavaanModel = lavaanModel, 
                                                penalty = penalty, 
                                                weights = weights_s, 
                                                tuningParameters = tp, 
                                                method = method,
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
        transformations = regularizedSEM_full@transformations,
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