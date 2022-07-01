#' cappedL1
#' 
#' This function allows for regularization of models built in lavaan with the
#' capped l1 penalty. This penalty is similar to the lasso but has an additional
#' tuning parameter (theta). All parameter estimates above this penalty will be 
#' penalized with constant theta. 
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
#' set.seed(123)
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
#' regsem <- cappedL1(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1),
#'   thetas = .1, # we "protect" paramter estimates which are above .1 in absolute
#'   # value. These parameters contribute to the -2log-Likelihood as usual, 
#'   # but are penalized as if they were of value .1
#'   control = controlIsta())
#' 
#' # use the plot-function to plot the regularized parameters:
#' plot(regsem)
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' 
#' # AIC and BIC:
#' AIC(regsem)
#' BIC(regsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(regsem, criterion = "AIC")
#' coef(regsem, criterion = "BIC")
#' 
#' # Note: theta is also a tuning parameter, so you may want to test
#' # multiple theta as well
#' regsem <- cappedL1(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1),
#'   thetas = seq(0.1,1,.1),
#'   control = controlIsta()) 
#' 
#' # The best parameters can also be extracted with:
#' coef(regsem, criterion = "AIC")
#' coef(regsem, criterion = "BIC")
#' 
#' # cross-validation
#' cv <- cv4cappedL1(regularizedSEM = regsem, 
#'                   k = 5)
#' # currently lacks a correct definition of coef, etc. Elements 
#' # can be accessed as follows:
#' cv@cvfits
#' @export
cappedL1 <- function(lavaanModel,
                     regularized,
                     lambdas,
                     thetas,
                     modifyModel = lessSEM::modifyModel(),
                     control = controlIsta()){
  if(any(thetas <= 0)) stop("Theta must be > 0")
  SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                 transformVariances = TRUE,
                                 whichPars = "est",
                                 addMeans = modifyModel$addMeans, 
                                 activeSet = modifyModel$activeSet,
                                 dataSet = modifyModel$dataSet)
  
  weights <- lessSEM::getParameters(SEM, raw = FALSE)
  weights[] <- 0
  weights[regularized] <- 1
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  result <- cappedElasticNet(
    lavaanModel = lavaanModel,
    weights = weights,
    lambdas = lambdas,
    alphas = 1,
    thetas = thetas,
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}

#' cappedElasticNet
#' 
#' This function allows for regularization of models built in lavaan with a combination
#' of ridge and capped l1 penalty. THIS IS EXPERIMENTAL AND SHOULD NOT BE USED!
#' This capped L1 penalty is similar to the lasso but has an additional
#' tuning parameter (theta; see ?cappedL1). All parameter estimates above this penalty will be 
#' penalized with constant theta. 
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
#' @param weights labeled vector with weights for each of the parameters in the 
#' model. If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' in [0,1]. 0 = ridge, 1 = lasso.
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta().
#' @md
#' @examples 
#' # see ?cappedL1, the use of cappedElasticNet is currently not recommended
#' @export
cappedElasticNet <- function(lavaanModel,
                             weights,
                             lambdas,
                             alphas,
                             thetas,
                             modifyModel = lessSEM::modifyModel(),
                             control = controlIsta()){
  
  inputArguments <- as.list(environment())
  
  if(!is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  
  if(!is(lavaanModel, "lavaan"))
    stop("lavaanModel must be of class lavaan")
  
  if(lavaanModel@Options$estimator != "ML") 
    stop("lavaanModel must be fit with ml estimator.")
  
  rawData <- try(lavaan::lavInspect(lavaanModel, "data"))
  if(is(rawData, "try-error")) 
    stop("Error while extracting raw data from lavaanModel. Please fit the model using the raw data set, not the covariance matrix.")
  N <- nrow(rawData)
  
  ### initialize model ####
  startingValues <- control$startingValues
  if(any(startingValues == "est")){
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "est",
                                   addMeans = modifyModel$addMeans, 
                                   activeSet = modifyModel$activeSet,
                                   dataSet = modifyModel$dataSet)
  }else if(any(startingValues == "start")){
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start",
                                   addMeans = modifyModel$addMeans, 
                                   activeSet = modifyModel$activeSet,
                                   dataSet = modifyModel$dataSet)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(lessSEM::getLavaanParameters(lavaanModel))))
      stop("Parameter names of startingValues do not match those of the lavaan object. See lessSEM::getLavaanParameters(lavaanModel).")
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                   transformVariances = TRUE,
                                   whichPars = "start", 
                                   fit = FALSE,
                                   addMeans = modifyModel$addMeans, 
                                   activeSet = modifyModel$activeSet,
                                   dataSet = modifyModel$dataSet)
    SEM <- lessSEM:::setParameters(SEM = SEM, labels = names(startingValues), value = startingValues, raw = FALSE)
    SEM <- try(lessSEM:::fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) 
      stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to elasticNet. See e.g., ?controlIsta for more information.")
  }
  
  # get parameters
  startingValues <- lessSEM:::getParameters(SEM, raw = TRUE)
  rawParameters <- lessSEM:::getParameters(SEM, raw = TRUE)
  
  # make sure that the weights are in the correct order
  if(is.null(names(weights))) stop("weights must have the same names as the parameters")
  if(length(weights) != length(rawParameters)) stop("weights must be of the same length as the parameter vector.")
  if(any(!is.numeric(weights))) stop("weights must be numeric")
  weights <- weights[names(rawParameters)]
  
  #### prepare regularized model object ####
  
  controlIntern <- list(
    L0 = control$L0,
    eta = control$eta,
    accelerate = control$accelerate,
    maxIterOut = control$maxIterOut,
    maxIterIn = control$maxIterIn,
    breakOuter = N*control$breakOuter,
    convCritInner = control$convCritInner,
    sigma = control$sigma,
    stepSizeInheritance = control$stepSizeInheritance,
    verbose = control$verbose
  )
  
  regularizedModel <- new(istaCappedL1, 
                          weights, 
                          controlIntern)
  
  #### define tuning parameters and prepare fit results ####
  
  tuningGrid <- expand.grid(
    "theta" = thetas,
    "lambda" = lambdas, 
    "alpha" = alphas)
  
  fits <- data.frame(
    tuningGrid,
    "m2LL" = NA,
    "regM2LL"= NA,
    "nonZeroParameters" = NA,
    "convergence" = NA
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,
                                             nrow = nrow(tuningGrid), 
                                             ncol = length(rawParameters)))
  colnames(parameterEstimates) <- names(rawParameters)
  parameterEstimates <- cbind(
    tuningGrid,
    parameterEstimates
  )
  
  Hessians <- list(NULL)
  
  #### print progress ####
  if(control$verbose == 0){
    progressbar = txtProgressBar(min = 0, 
                                 max = nrow(tuningGrid), 
                                 initial = 0, 
                                 style = 3)
  }
  
  #### Iterate over all tuning parameter combinations and fit models ####
  
  for(it in 1:nrow(tuningGrid)){
    if(control$verbose == 0){
      setTxtProgressBar(progressbar,it)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tuningGrid),"]\n"))
    }
    
    lambda <- tuningGrid$lambda[it]
    alpha <- tuningGrid$alpha[it]
    theta <- tuningGrid$theta[it]
    
    result <- try(regularizedModel$optimize(rawParameters,
                                            SEM,
                                            theta,
                                            lambda,
                                            alpha)
    )
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(rawParameters[weights[names(rawParameters)] != 0] == 0)
    fits$regM2LL[it] <- result$fit
    fits$convergence[it] <- result$convergence
    
    SEM <- lessSEM::setParameters(SEM, 
                                  names(rawParameters), 
                                  values = rawParameters, 
                                  raw = TRUE)
    fits$m2LL[it] <- SEM$fit()
    transformedParameters <- lessSEM:::getParameters(SEM,
                                                     raw = FALSE)
    parameterEstimates[it, names(rawParameters)] <- transformedParameters[names(rawParameters)]
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning("Fit for lambda = ",
              lambda, "alpha = ", 
              alpha, "theta = ", theta,
              " resulted in Error!")
      
      SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                     transformVariances = TRUE,
                                     whichPars = startingValues,
                                     addMeans = control$addMeans)
      
      
    }
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("regularizedSEM",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}

#' cv4cappedL1
#' 
#' Exact cross-validation for models of class regularizedSEM which have been fitted with cappedL1 penalty
#' These models can be fit with cappedL1() (see ?cappedL1)
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
#' # see ?cappedL1
#' 
#' @export
cv4cappedL1 <- function(regularizedSEM, 
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
  return(cv4cappedElasticNet(regularizedSEM, 
                             k,
                             dataSet = dataSet,
                             scalingFunction = scalingFunction,
                             scalingArguments = scalingArguments,
                             reweigh = FALSE,
                             returnSubsetParameters = FALSE))
  
}



#' cv4cappedElasticNet
#' 
#' Exact cross-validation for models of class regularizedSEM which have been fitted a combination
#' of ridge and capped l1 penalty. THIS IS EXPERIMENTAL AND SHOULD NOT BE USED!
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
#' # see ?cv4cappedL1
#' @export
cv4cappedElasticNet <- function(regularizedSEM, 
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
  
  tuningParameters <- data.frame(theta = fits$theta,
                                 lambda = fits$lambda, 
                                 alpha = fits$alpha
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
          method = "ista", 
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
    
    regularizedSEM_s <- lessSEM::cappedElasticNet(
      lavaanModel = regularizedSEM@inputArguments$lavaanModel,
      weights = weights_s,
      lambdas = regularizedSEM@inputArguments$lambdas, 
      alphas = regularizedSEM@inputArguments$alphas, 
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