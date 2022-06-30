#### This script provides functions to perform exact cross-validation for different types of models ####

#' cv4lasso
#' 
#' Exact cross-validation for models of class regularizedSEM which have been fitted with lasso penalty
#' These models can be fit with lasso() (see ?lasso)
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
#' library(lessSEM)
#' 
#' @export
cv4lasso <- function(regularizedSEM, 
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
  return(cv4elasticNet(regularizedSEM, 
                       k,
                       dataSet = dataSet,
                       scalingFunction = scalingFunction,
                       scalingArguments = scalingArguments,
                       reweigh = FALSE,
                       returnSubsetParameters = FALSE))
  
}

#' cv4adaptiveLasso
#' 
#' Exact cross-validation for models of class regularizedSEM which have been 
#' fitted with adaptive lasso penalty
#' These models can be fit with adaptiveLasso() (see ?adaptiveLasso)
#' in this package. A difficulty when cross-validating the adaptive lasso is
#' that the weights are often based on the maximum likelihood estimates of 
#' the full sample. Using the same weights for the subsets of the cross-validation
#' can undermine the independence of these subsets. We therefore recommend reweighing
#' the subsets. If you used the default weights in adaptiveLasso (setting weights = NULL)
#' you can set reweigh = TRUE and cv4adaptiveLasso will automatically re-compute the
#' MLE for each training set. These MLE_training_set will then be used to construct the weights
#' as weight = 1/abs(MLE_training_set). For unregularized parameters, weights will be set to zero.
#' You can also pass custom weights. To this end, you can pass a matrix with k rows and as many
#' columns as there are parameters in your model to the reweigh argument. Note that you 
#' should then also pass specific data sets instead of using e.g. k = 5.
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
#' @param reweigh this is used when weight are not all 1 or 0 (in adaptive lasso). See ?cv4adaptiveLasso for details
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @examples
#' library(lessSEM)
#' 
#' @export
cv4adaptiveLasso <- function(regularizedSEM, 
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
  return(cv4elasticNet(regularizedSEM, 
                       k,
                       dataSet = dataSet,
                       scalingFunction = scalingFunction,
                       scalingArguments = scalingArguments,
                       reweigh = reweigh,
                       returnSubsetParameters = FALSE))
  
}

#' cv4ridge
#' 
#' Exact cross-validation for models of class regularizedSEM which have been fitted with ridge penalty
#' These models can be fit with ridge() (see ?ridge)
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
#' library(lessSEM)
#' 
#' @export
cv4ridge <- function(regularizedSEM, 
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
  return(cv4elasticNet(regularizedSEM, 
                       k,
                       dataSet = dataSet,
                       scalingFunction = scalingFunction,
                       scalingArguments = scalingArguments,
                       reweigh = FALSE,
                       returnSubsetParameters = FALSE))
  
}

#' cv4elasticNet
#' 
#' Exact cross-validation for models of class regularizedSEM which have been fitted with elastic net penalty
#' These models can be fit with elasticNet() (see ?elasticNet)
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
#' @param reweigh this is used when weight are not all 1 or 0 (in adaptive lasso). See ?cv4adaptiveLasso for details
#' @param returnSubsetParameters if set to TRUE, the parameter estimates of the individual cross-validation training sets will be returned
#' @examples
#' library(lessSEM)
#' 
#' # Let's first set up a regularized model. The following steps are
#' # explained in detail in ?lessSEM::regularizeSEM
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
#' ## The exact cross-validation can be computed with:
#' CV <- aCV4regularizedSEM(regularizedSEM = regsem,
#'                           k = 5
#' )
#' # let's plot the parameter values and the corresponding leave-one-out fit:
#' plot(CV)
#' 
#' # the best parameters can be extracted with
#' coef(CV)
#' 
#' # we can also use the one standard deviation rule:
#' coef(CV, rule = "1sd")
#' 
#' # To see which person ended up in which sample, use:
#' CV@subsets
#' 
#' #### Advanced ####
#' # If you are interested in the parameter estimates of each sub-sample,
#' # you must re-run the computation with returnSubsetParameters set to TRUE.
#' # This is disabled by default because it may take a lot of disk space
#' CV <- CV4regularizedSEM(regularizedSEM = regsem,
#'                           k = 5,
#'                           returnSubsetParameters = TRUE
#' )
#' 
#' # the parameter are returned in a 3D-array. The subsets are in the
#' # rows, the parameters in the columns and the tuning parameters in
#' # the third dimension. To access the elements, use:
#' # Access elements for:
#' dimnames(CV@subsetParameters)[[3]][1]
#' CV@subsetParameters[,,1] # first lambda value
#' # Access elements for:
#' dimnames(CV@subsetParameters)[[3]][2]
#' CV@subsetParameters[,,2] # second lambda value
#' 
#' ## Instead of letting lessSEM create the subsets, you can
#' # also pass your own subsets. As an example:
#' subsets <- lessSEM:::createSubsets(N = nrow(dataset), k  = 10)
#' 
#' CV <- CV4regularizedSEM(regularizedSEM = regsem,
#'                           k = subsets
#' )
#' CV@cvfitsDetails
#' @export
cv4elasticNet <- function(regularizedSEM, 
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
  if(any(!weights %in% c(0,1))) message("Automatic cross-validation for adaptiveLasso requested. Note that using weights which are based on the full sample may undermine cross-validation. Use reweigh = TRUE to re-weight each subset with the inverse of the absolute MLE. Alternatively, pass a matrix as argument reweigh with weights for each subset.")
  
  tuningParameters <- data.frame(lambda = fits$lambda, 
                                 alpha = fits$alpha)
  
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
    if(any(!weights %in% c(0,1)) && any(reweigh)){
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
        param_s <- coef(MLE_s, lambda = 0, alpha = 1)
        weights_s <- 1/abs(param_s)
        # set unregularized to 0:
        weights_s[regularizedSEM@inputArguments$weights == 0] <- 0
      }
      
    }else{
      weights_s <- regularizedSEM@inputArguments$weights
    }
    
    regularizedSEM_s <- lessSEM::elasticNet(
      lavaanModel = regularizedSEM@inputArguments$lavaanModel,
      weights = weights_s,
      lambdas = regularizedSEM@inputArguments$lambdas, 
      nLambdas = NULL, 
      alphas = regularizedSEM@inputArguments$alphas, 
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
