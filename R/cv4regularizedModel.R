#' cv4regularizedSEM
#' 
#' cross-validation for regularized structural equation models
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
#' @param scalingArguments the second argument passed to scalingFunction.
#' @param reweigh this is used for the adaptive lasso. When the weights are based on the full
#' sample, this may undermine the cross-validation. Set reweigh = TRUE to create new
#' weigths for each subset. This will use the default weights (inverse of MLE). Alternatively,
#' you can pass a matrix with k rows and nParameters columns with weights.
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' 
#' @examples 
#' library(lavaan)
#' library(lessSEM)
#' set.seed(123)
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
#' regsem <- lasso(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1),
#'   control = controlIsta())
#' )
#' 
#' # standard cross-validation
#' cv <- try(cv4regularizedSEM(regularizedSEM = regsem,
#'                               k = 5)
#' )
#' 
#' # IMPORTANT: If your data is scaled prior to the cross-validation, 
#' # there are dependencies between your train and test set. To avoid this,
#' # you can use the scale arguments to standardize the training and testing sets
#' # separately.
#' # First, set scaleData = TRUE
#' # Second, provide a scalingFunction. This function MUST take two arguments:
#' # the data set and the scalingArguments.
#' # scalingArguments can be used to pass anything you need to the scalingFunction.
#' # Example:
#' scalingFunction <-  function(x, scalingArguments) {
#'   scale(x, 
#'         center = scalingArguments$center, 
#'         scale = scalingArguments$scale)
#' }
#' scalingArguments = list("scale" = TRUE, 
#'                         "center" = TRUE)
#' 
#' cv <- cv4regularizedSEM(regularizedSEM = regsem,
#'                           k = 5,
#'                           scaleData = TRUE, 
#'                           scalingFunction = scalingFunction, 
#'                           scalingArguments = scalingArguments )
#' 
#' # In adaptive lasso regularization, using adaptive lasso weights which 
#' # are based on the full sample may also result in false cross-validation.
#' # Here, you can use the reweigh argument to re-compute the weights for each
#' # train set. Note that these will always be the inverse of the absolute values
#' # of the maximum likelihood estimates. Alternatively, you can also pass a matrix
#' # with k rows and nParameter columns with weights
#' 
#' regsem <-  adaptiveLasso(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1),
#'   control = controlIsta()
#' ) 
#' 
#' cv <- cv4regularizedSEM(regularizedSEM = regsem,
#'                           k = 5, 
#'                           reweigh = TRUE)
#' @export
cv4regularizedSEM <- function(regularizedSEM, 
                                k, 
                                dataSet = NULL,
                                scaleData = FALSE,
                                scalingFunction = function(dataSet,scalingArguments) 
                                  scale(x = dataSet, 
                                        center = scalingArguments$center, 
                                        scale = scalingArguments$scale),
                                scalingArguments = list("center" = TRUE, 
                                                        "scale" = TRUE),
                                reweigh = FALSE,
                                returnSubsetParameters = FALSE){
  
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  
  lavaanData <- try(lavaan::lavInspect(regularizedSEM@inputArguments$lavaanModel,
                                       "data"))
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
  
  if(regularizedSEM@penalty == "adaptiveLasso") 
    message("Automatic cross-validation for adaptiveLasso requested. Note that using weights which are based on the full sample may undermine cross-validation. Use reweigh = TRUE to re-weight each subset with the inverse of the absolute MLE. Alternatively, pass a matrix as argument reweigh with weights for each subset.")
  
  tuningParameters <- parameters[,!colnames(parameters)%in%regularizedSEM@parameterLabels, drop = FALSE]
  
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
    
    # check weights for adaptive Lasso
    if(regularizedSEM@penalty == "adaptiveLasso" && any(reweigh)){
      if(is.matrix(reweigh)){
        weights_s <- reweigh[s,]
      }else{
        # use default;
        print(paste0("Computing new weights for sample ", s, "."))
        
        # optimize model: We set lambda = 0, so we get the MLE
        MLE_s <- lessSEM::lasso(
          lavaanModel = regularizedSEM@inputArguments$lavaanModel,
          lambdas = 0, 
          regularized = regularizedSEM@regularized,
          method = regularizedSEM@inputArguments$method, 
          modifyModel = modifyModel,
          control = control_s
        )
        # MLE:
        param_s <- unlist(MLE_s@parameters[,MLE_s@parameterLabels])
        weights_s <- 1/abs(param_s)
        # set unregularized to 0:
        if(!is.numeric(regularizedSEM@inputArguments$weights)){
          weights_s[!names(weights_s) %in% regularizedSEM@inputArguments$weights] <- 0
        }else{
          weights_s[regularizedSEM@inputArguments$weights == 0] <- 0
        }
      }
      
    }else{
      weights_s <- regularizedSEM@inputArguments$weights
    }
    
    regularizedSEM_s <- regularizeSEMInternal(lavaanModel = regularizedSEM@inputArguments$lavaanModel, 
                                              penalty = regularizedSEM@penalty, 
                                              weights = weights_s, 
                                              tuningParameters = tuningParameters, 
                                              method = regularizedSEM@inputArguments$method,
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


#' cv4regularizedSEMApprox
#' 
#' cross-validation for regularized structural equation models with approximated penalty function
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
#' @param scalingArguments the second argument passed to scalingFunction.
#' @param reweigh this is used for the adaptive lasso. When the weights are based on the full
#' sample, this may undermine the cross-validation. Set reweigh = TRUE to create new
#' weigths for each subset. This will use the default weights (inverse of MLE). Alternatively,
#' you can pass a matrix with k rows and nParameters columns with weights.
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @export
cv4regularizedSEMApprox <- function(regularizedSEM, 
                                    k, 
                                    dataSet = NULL,
                                    scaleData = FALSE,
                                    scalingFunction = function(dataSet,scalingArguments) 
                                      scale(x = dataSet, 
                                            center = scalingArguments$center, 
                                            scale = scalingArguments$scale),
                                    scalingArguments = list("center" = TRUE, 
                                                            "scale" = TRUE),
                                    reweigh = FALSE,
                                    returnSubsetParameters = FALSE){
  
  if(!is(regularizedSEM, "regularizedSEM")){
    stop("regularizedSEM must be of class regularizedSEM")
  }
  
  lavaanData <- try(lavaan::lavInspect(regularizedSEM@inputArguments$lavaanModel,
                                       "data"))
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
  epsilon <- regularizedSEM@inputArguments$epsilon
  tau <- regularizedSEM@inputArguments$tau
  
  if(any(!fits$alpha %in% c(0,1))) 
    message("Automatic cross-validation for adaptiveLasso requested. Note that using weights which are based on the full sample may undermine cross-validation. Use reweigh = TRUE to re-weight each subset with the inverse of the absolute MLE. Alternatively, pass a matrix as argument reweigh with weights for each subset.")
  
  tuningParameters <- parameters[,!colnames(parameters)%in%regularizedSEM@parameterLabels, drop = FALSE]
  
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
    
    # check weights for adaptive Lasso
    if(any(!fits$alpha %in% c(0,1)) && any(reweigh)){
      if(is.matrix(reweigh)){
        weights_s <- reweigh[s,]
      }else{
        # use default;
        print(paste0("Computing new weights for sample ", s, "."))
        
        # optimize model: We set lambda = 0, so we get the MLE
        MLE_s <- lessSEM::smoothLasso(
          lavaanModel = regularizedSEM@inputArguments$lavaanModel,
          lambdas = 0,
          regularized = regularizedSEM@regularized,
          tau = tau,
          epsilon = epsilon,
          modifyModel = modifyModel,
          control = control_s
        )
        # MLE:
        param_s <- unlist(MLE_s@parameters[,MLE_s@parameterLabels])
        weights_s <- 1/abs(param_s)
        # set unregularized to 0:
        if(!is.numeric(regularizedSEM@inputArguments$weights)){
          weights_s[!names(weights_s) %in% regularizedSEM@inputArguments$weights] <- 0
        }else{
          weights_s[regularizedSEM@inputArguments$weights == 0] <- 0
        }
      }
      
    }else{
      weights_s <- regularizedSEM@inputArguments$weights
    }
    
    regularizedSEM_s <- smoothElasticNet(lavaanModel = regularizedSEM@inputArguments$lavaanModel, 
                                         weights = weights_s, 
                                         lambdas = unique(tuningParameters$lambda), 
                                         alphas = unique(tuningParameters$alpha),
                                         modifyModel = modifyModel,
                                         tau = tau,
                                         epsilon = epsilon,
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