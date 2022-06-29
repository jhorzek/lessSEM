#' smoothLasso
#' 
#' This function allows for regularization of models built in lavaan with the
#' smoothed lasso penalty. The returned object is an S4 class; its elements can be accessed
#' with the "@" operator (see examples). We don't recommend using this function.
#' Use lasso() instead.
#' 
#' For more details, see:
#' 
#' 1. Lee, S.-I., Lee, H., Abbeel, P., & Ng, A. Y. (2006). 
#' Efficient L1 Regularized Logistic Regression. Proceedings of the 
#' Twenty-First National Conference on Artificial Intelligence (AAAI-06), 401–408.
#' 2. Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). 
#'    Regularized Structural Equation Modeling. Structural Equation Modeling: 
#'    A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793    
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, linr can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param tau parameters below threshold tau will be seen as zeroed
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @md
#' @examples 
#' library(linr)
#' 
#' # Identical to regsem, linr builds on the lavaan
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
#' regsem <- smoothLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   epsilon = 1e-10,
#'   tau = 1e-4,
#'   lambdas = seq(0,1,length.out = 50))
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
#' @export
smoothLasso <- function(lavaanModel,
                        regularized,
                        lambdas = NULL,
                        nLambdas = NULL,
                        epsilon,
                        tau,
                        control = controlBFGS()){
  
  weights <- getLavaanParameters(lavaanModel)
  weights[] <- 0
  weights[regularized] <- 1
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  result <- smoothElasticNet(
    lavaanModel = lavaanModel,
    weights = weights,
    lambdas = lambdas,
    nLambdas = nLambdas,
    alphas = 1,
    epsilon = epsilon, 
    tau = tau,
    control = control
  )
  return(result)
  
}

#' smoothAdaptiveLasso
#' 
#' This function allows for regularization of models built in lavaan with the
#' smooth adaptive lasso penalty. The returned object is an S4 class; its elements can be accessed
#' with the "@" operator (see examples).
#' 
#' For more details, see:
#' 
#' 1. Zou, H. (2006). The Adaptive Lasso and Its Oracle Properties.
#'    Journal of the American Statistical Association, 101(476), 1418–1429.
#'    https://doi.org/10.1198/016214506000000735
#' 2. Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). 
#'    Regularized Structural Equation Modeling. Structural Equation Modeling: 
#'    A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 3. Lee, S.-I., Lee, H., Abbeel, P., & Ng, A. Y. (2006). 
#' Efficient L1 Regularized Logistic Regression. Proceedings of the 
#' Twenty-First National Conference on Artificial Intelligence (AAAI-06), 401–408.
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param weights labeled vector with weights for each of the parameters in the 
#' model. If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object. If set to NULL,
#' the default weights will be used: the inverse of the absolute values of
#' the unregularized parameter estimates
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, linr can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param tau parameters below threshold tau will be seen as zeroed
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @md
#' @examples 
#' library(linr)
#' 
#' # Identical to regsem, linr builds on the lavaan
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
#' # names of the regularized parameters:
#' regularized = paste0("l", 6:15)
#' 
#' # define adaptive lasso weights:
#' # We use the inverse of the absolute unregularized parameters
#' # (this is the default in adaptiveLasso and can also specified
#' # by setting weights = NULL)
#' weights <- 1/abs(getLavaanParameters(lavaanModel))
#' weights[!names(weights) %in% regularized] <- 0
#' 
#' regsem <- smoothAdaptiveLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   regularized = regularized,
#'   weights = weights,
#'   epsilon = 1e-10,
#'   tau = 1e-4,
#'   lambdas = seq(0,1,length.out = 50))
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
#' @export
smoothAdaptiveLasso <- function(lavaanModel,
                          regularized,
                          weights = NULL,
                          lambdas = NULL,
                          nLambdas = NULL,
                          epsilon,
                          tau,
                          control = controlBFGS()){
  if(is.null(weights)){
    weights <- 1/abs(getLavaanParameters(lavaanModel))
    weights[!names(weights) %in% regularized] <- 0
  }
  
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  result <- smoothElasticNet(
    lavaanModel = lavaanModel,
    weights = weights,
    lambdas = lambdas,
    nLambdas = nLambdas,
    alphas = 1,
    epsilon = epsilon, 
    tau = tau,
    control = control
  )
  return(result)
  
}

#' ridgeBfgs
#' 
#' This function allows for regularization of models built in lavaan with the
#' ridge penalty. Its elements can be accessed
#' with the "@" operator (see examples).
#' 
#' For more details, see:
#' 
#' 1. Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). 
#'    Regularized Structural Equation Modeling. Structural Equation Modeling: 
#'    A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 2. Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood 
#'    Method for Structural Equation Modeling. Psychometrika, 82(2),
#'    329–354. https://doi.org/10.1007/s11336-017-9566-9
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlBFGS for more details.
#' @md
#' @examples 
#' library(linr)
#' 
#' # Identical to regsem, linr builds on the lavaan
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
#' # names of the regularized parameters:
#' regularized = paste0("l", 6:15)
#' 
#' regsem <- ridgeBfgs(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   regularized = regularized,
#'   lambdas = seq(0,1,length.out = 50))
#' 
#' plot(regsem)
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' @export
ridgeBfgs <- function(lavaanModel,
                  regularized,
                  lambdas = NULL,
                  control = controlBFGS()){
  
  weights <- getLavaanParameters(lavaanModel)
  weights[] <- 0
  weights[regularized] <- 1
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  result <- smoothElasticNet(
    lavaanModel = lavaanModel,
    weights = weights,
    lambdas = lambdas,
    nLambdas = NULL,
    alphas = 0,
    epsilon = 0, # ridge is already smooth
    tau = 0,
    control = control
  )
  return(result)
  
}

#' smoothElasticNet
#' 
#' This function allows for regularization of models built in lavaan with the
#' smooth elastic net penalty. Its elements can be accessed
#' with the "@" operator (see examples).
#' 
#' For more details, see:
#' 
#' 1. Zou, H., & Hastie, T. (2005). Regularization and variable selection via 
#' the elastic net. Journal of the Royal Statistical Society: 
#' Series B, 67(2), 301–320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
#' for the details of this regularization technique.
#' 2. Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). 
#' Regularized Structural Equation Modeling. Structural Equation Modeling: 
#' A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 3. Lee, S.-I., Lee, H., Abbeel, P., & Ng, A. Y. (2006). 
#' Efficient L1 Regularized Logistic Regression. Proceedings of the 
#' Twenty-First National Conference on Artificial Intelligence (AAAI-06), 401–408.
#' 
#' @param lavaanModel model of class lavaan 
#' @param weights labeled vector with weights for each of the parameters in the 
#' model. If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, linr can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' in [0,1]. 0 = ridge, 1 = lasso.
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param tau parameters below threshold tau will be seen as zeroed
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlBFGS for more details.
#' @md
#' @examples 
#' library(linr)
#' 
#' # Identical to regsem, linr builds on the lavaan
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
#' # names of the regularized parameters:
#' regularized = paste0("l", 6:15)
#' 
#' # define adaptive lasso weights:
#' # We use the inverse of the absolute unregularized parameters
#' # (this is the default in adaptiveLasso and can also specified
#' # by setting weights = NULL)
#' weights <- getLavaanParameters(lavaanModel)
#' weights[] <- 0
#' weights[regularized] <- 1
#' 
#' regsem <- smoothElasticNet(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   weights = weights,
#'   epsilon = 1e-10,
#'   tau = 1e-4,
#'   lambdas = seq(0,1,length.out = 50),
#'   alphas = seq(0,1,length.out = 4))
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' @export
smoothElasticNet <- function(lavaanModel,
                       weights,
                       lambdas = NULL,
                       nLambdas = NULL,
                       alphas,
                       epsilon, 
                       tau,
                       control = controlBFGS()){
  
  inputArguments <- as.list(environment())
  
  if(any(alphas != 0)) warning("Optimization of non-smooth functions with bfgs is not recommended.")
  
  if(is.null(lambdas) && is.null(nLambdas)) 
    stop("Specify lambdas or nLambdas")
  if(!is.null(lambdas) && !is.null(nLambdas))
    stop("Specify either lambdas or nLambdas, not both")
  if(!is.null(nLambdas) && any(alphas != 1))
    stop("nLambdas only supported for lasso type penalty (alpha = 1).")
  
  if(!is(control, "controlBFGS")) 
    stop("control must be of class controlBFGS See ?controlBFGS")
  
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
    SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "est",
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
  }else if(any(startingValues == "start")){
    SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "start",
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(linr::getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See linr::getLavaanParameters(lavaanModel).")
    SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "start", 
                                fit = FALSE,
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
    SEM <- linr:::setParameters(SEM = SEM, labels = names(startingValues), value = startingValues, raw = FALSE)
    SEM <- try(linr:::fit(SEM))
    if(is(SEM, "try-error") || !is.finite(SEM$m2LL)) 
      stop("Infeasible starting values.")
    
  }else{
    stop("Invalid startingValues passed to elasticNet. See e.g., ?controlIsta for more information.")
  }
  
  # get parameters
  startingValues <- linr:::getParameters(SEM, raw = TRUE)
  rawParameters <- linr:::getParameters(SEM, raw = TRUE)
  
  # make sure that the weights are in the correct order
  if(is.null(names(weights))) stop("weights must have the same names as the parameters")
  if(length(weights) != length(rawParameters)) stop("weights must be of the same length as the parameter vector.")
  if(any(!is.numeric(weights))) stop("weights must be numeric")
  weights <- weights[names(rawParameters)]
  
  #### bfgs requires an initial Hessian ####
  initialHessian <- control$initialHessian
  if(is.matrix(initialHessian) && nrow(initialHessian) == length(rawParameters) && ncol(initialHessian) == length(rawParameters)){
    
    if(!all(rownames(initialHessian) %in% names(linr::getLavaanParameters(lavaanModel))) ||
       !all(colnames(initialHessian) %in% names(linr::getLavaanParameters(lavaanModel)))
    ) stop("initialHessian must have the parameter names as rownames and colnames. See linr::getLavaanParameters(lavaanModel).")
    
  }else if(any(initialHessian == "compute")){
    
    initialHessian <- linr:::getHessian(SEM = SEM, raw = TRUE)
    
  }else if(any(initialHessian == "scoreBased")){
    
    scores <- linr:::getScores(SEM = SEM, raw = TRUE)
    FisherInformation <- matrix(0, nrow = ncol(scores), ncol = ncol(scores))
    rownames(FisherInformation) <- colnames(FisherInformation) <- colnames(scores)
    for(score in 1:nrow(scores)){
      FisherInformation <- FisherInformation + t(-.5*scores[score,, drop = FALSE]) %*%(-.5*scores[score,, drop = FALSE]) # we are using the -2 log-Likelihood
    }
    
    initialHessian <- -2*(-FisherInformation) # negative log-likelihood
    # make symmetric; just in case...
    initialHessian <- .5*(initialHessian + t(initialHessian))
    while(any(eigen(initialHessian, only.values = TRUE)$values < 0)){
      diag(initialHessian) <- diag(initialHessian) + 1e-4
    }
    
  }else if(length(initialHessian) == 1 && is.numeric(initialHessian)){
    initialHessian <- diag(initialHessian,length(rawParameters))
    rownames(initialHessian) <- names(rawParameters)
    colnames(initialHessian) <- names(rawParameters)
  }else{
    stop("Invalid initialHessian passed to BFGS. See ?controlBFGS for more information.")
  }
  
  control$initialHessian <- initialHessian
  
  #### prepare regularized model object ####
  
  controlIntern <- list(
    epsilon = epsilon,
    initialHessian = control$initialHessian,
    stepSize = control$stepSize,
    sigma = control$sigma,
    gamma = control$gamma,
    maxIterOut = control$maxIterOut,
    maxIterIn = control$maxIterIn,
    maxIterLine = control$maxIterLine,
    breakOuter = N*control$breakOuter,
    breakInner = N*control$breakInner,
    convergenceCriterion = control$convergenceCriterion, 
    verbose = control$verbose
  )
  
  regularizedModel <- new(bfgsEnet, 
                          weights, 
                          controlIntern)
  
  #### define tuning parameters and prepare fit results ####
  ## get max lambda ##
  if(!is.null(nLambdas)){
    message(paste0(
      "Automatically selecting the maximal lambda value.\n",
      "Note: This may fail if a model with all regularized parameters set to zero is not identified.")
    )
    
    maxLambda <- getMaxLambda_C(regularizedModel = regularizedModel,
                                SEM = SEM,
                                rawParameters = rawParameters,
                                weights = weights,
                                N = N)
    lambdas <- seq(0,
                   maxLambda,
                   length.out = nLambdas)
    
  }
  
  tuningGrid <- expand.grid("lambda" = lambdas, 
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
  
  if(control$saveHessian){
    Hessians <- list(
      "lambda" = tuningGrid$lambda,
      "alpha" = tuningGrid$alpha,
      "Hessian" = lapply(1:nrow(tuningGrid), 
                         matrix, 
                         data= NA, 
                         nrow=nrow(initialHessian), 
                         ncol=ncol(initialHessian))
    )
  }else{
    Hessians <- list(NULL)
  }
  
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
    
    result <- try(regularizedModel$optimize(rawParameters,
                                            SEM,
                                            lambda,
                                            alpha)
    )
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(abs(rawParameters[weights[names(rawParameters)] != 0]) <= tau)
    fits$regM2LL[it] <- result$fit
    fits$convergence[it] <- result$convergence
    
    SEM <- linr::setParameters(SEM, 
                               names(rawParameters), 
                               values = rawParameters, 
                               raw = TRUE)
    fits$m2LL[it] <- SEM$fit()
    transformedParameters <- linr:::getParameters(SEM,
                                                  raw = FALSE)
    parameterEstimates[it, names(rawParameters)] <- transformedParameters[names(rawParameters)]
    
      Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning("Fit for lambda = ",
              lambda, "alpha = ", 
              alpha,
              " resulted in Error!")
      
      SEM <- linr:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                  transformVariances = TRUE,
                                  whichPars = startingValues,
                                  addMeans = control$addMeans)
      
      regularizedModel$setHessian(controlIntern$initialHessian)
      
    }else{
      
      if(control$saveHessian) Hessians$Hessian[[it]] <- result$Hessian
      
      # set Hessian for next iteration
      regularizedModel$setHessian(result$Hessian)
      
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