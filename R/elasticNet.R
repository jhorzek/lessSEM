#' lasso
#' 
#' This function allows for regularization of models built in lavaan with the
#' lasso penalty. The returned object is an S4 class; its elements can be accessed
#' with the "@" operator (see examples).
#' 
#' For more details, see:
#' 
#' 1. Tibshirani, R. (1996). Regression Shrinkage and Selection via the Lasso. 
#'    Journal of the Royal Statistical Society. Series B (Methodological), 58(1), 267–288.
#' 2. Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). 
#'    Regularized Structural Equation Modeling. Structural Equation Modeling: 
#'    A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793    
#' 3. Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood 
#'    Method for Structural Equation Modeling. Psychometrika, 82(2),
#'    329–354. https://doi.org/10.1007/s11336-017-9566-9
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
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
#' regsem <- lasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   # in case of lasso and adaptive lasso, we can specify the number of lambda
#'   # values to use. lessSEM will automatically find lambda_max and fit
#'   # models for nLambda values between 0 and lambda_max. For the other
#'   # penalty functions, lambdas must be specified explicitly
#'   nLambdas = 50)
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
#' #### Advanced ###
#' # Switching the optimizer # 
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' regsemGlmnet <- lasso(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   nLambdas = 50,
#'   method = "glmnet",
#'   control = controlGlmnet())
#' 
#' # Note: The results are basically identical:
#' regsemGlmnet@parameters - regsem@parameters
#' 
#' ## The fitted model can then be used as basis for an approximate cross-validation
#' # (see ?lessSEM::acv4lasso) or approximate influence functions
#' # (see ?lessSEM::ai4lasso)
#' @export
lasso <- function(lavaanModel,
                  regularized,
                  lambdas = NULL,
                  nLambdas = NULL,
                  method = "ista", 
                  control = controlIsta()){
  SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                 transformVariances = TRUE,
                                 whichPars = "est",
                                 addMeans = control$addMeans, 
                                 activeSet = control$activeSet)
  
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
  
  result <- elasticNet(
    lavaanModel = lavaanModel,
    weights = weights,
    lambdas = lambdas,
    nLambdas = nLambdas,
    alphas = 1,
    method = method, 
    control = control
  )
  return(result)
  
}

#' adaptiveLasso
#' 
#' This function allows for regularization of models built in lavaan with the
#' adaptive lasso penalty. The returned object is an S4 class; its elements can be accessed
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
#' 3. Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood 
#'    Method for Structural Equation Modeling. Psychometrika, 82(2),
#'    329–354. https://doi.org/10.1007/s11336-017-9566-9
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
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
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
#' regsem <- adaptiveLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = regularized, 
#'   # set adaptive lasso weights
#'   weights = weights,
#'   # in case of lasso and adaptive lasso, we can specify the number of lambda
#'   # values to use. lessSEM will automatically find lambda_max and fit
#'   # models for nLambda values between 0 and lambda_max. For the other
#'   # penalty functions, lambdas must be specified explicitly
#'   nLambdas = 50)
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
#' #### Advanced ###
#' # Switching the optimizer # 
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' regsemGlmnet <- adaptiveLasso(
#'   lavaanModel = lavaanModel,
#'   regularized = regularized,
#'   weights = weights,
#'   nLambdas = 50,
#'   method = "glmnet",
#'   control = controlGlmnet())
#' 
#' # Note: The results are basically identical:
#' regsemGlmnet@parameters - regsem@parameters
#' @export
adaptiveLasso <- function(lavaanModel,
                  regularized,
                  weights = NULL,
                  lambdas = NULL,
                  nLambdas = NULL,
                  method = "ista", 
                  control = controlIsta()){
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
  
  result <- elasticNet(
    lavaanModel = lavaanModel,
    weights = weights,
    lambdas = lambdas,
    nLambdas = nLambdas,
    alphas = 1,
    method = method, 
    control = control
  )
  return(result)
  
}

#' ridge
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
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
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
#' regsem <- ridge(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   # for ridge regularization, the lambdas have to be defined
#'   # explicitly:
#'   lambdas = seq(0,1,length.out = 100))
#' 
#' # use the plot-function to plot the regularized parameters:
#' plot(regsem)
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' 
#' # Model selection should be done using cross-validation
#' cv <- cv4ridge(regularizedSEM = regsem, k = 5)
#' plot(cv)
#' coef(cv) # get best fitting model parameters
#' 
#' #### Advanced ###
#' # Switching the optimizer # 
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' regsemGlmnet <- ridge(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 100),
#'   method = "glmnet",
#'   control = controlGlmnet())
#' 
#' # Note: The results are basically identical:
#' regsemGlmnet@parameters - regsem@parameters
#' 
#' ## The fitted model can then be used as basis for an approximate cross-validation
#' # (see ?lessSEM::acv4ridge) or approximate influence functions
#' # (see ?lessSEM::ai4ridge)
#' @export
ridge <- function(lavaanModel,
                  regularized,
                  lambdas = NULL,
                  method = "ista", 
                  control = controlIsta()){
  
  SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                 transformVariances = TRUE,
                                 whichPars = "est",
                                 addMeans = control$addMeans, 
                                 activeSet = control$activeSet)
  
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
  
  result <- elasticNet(
    lavaanModel = lavaanModel,
    weights = weights,
    lambdas = lambdas,
    nLambdas = NULL,
    alphas = 0,
    method = method, 
    control = control
  )
  return(result)
  
}

#' elasticNet
#' 
#' This function allows for regularization of models built in lavaan with the
#' elastic net penalty. Its elements can be accessed
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
#' 3. Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood 
#' Method for Structural Equation Modeling. Psychometrika, 82(2),
#'  329–354. https://doi.org/10.1007/s11336-017-9566-9
#' 
#' @param lavaanModel model of class lavaan 
#' @param weights labeled vector with weights for each of the parameters in the 
#' model. If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' in [0,1]. 0 = ridge, 1 = lasso.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta() and controlGlmnet() functions.
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
#' # The implementation of elastic net is relatively flexible and also
#' # encompasses lasso, adaptive lasso, and ridge regularization.
#' # The function call is therefore slightly more complicated. If you 
#' # are only looking for lasso, adaptive lasso or ridge, there are 
#' # dedicated functions for that (see ?lasso, ?adaptiveLasso, ?ridge)
#' 
#' # elasticNet expects weights for each parameter. This way, we can
#' # decide which parameters are regularized (all parameters with weight = 0
#' # remain unregularized)
#' # names of the regularized parameters:
#' regularized = paste0("l", 6:15)
#' weights <- getLavaanParameters(lavaanModel)
#' weights[!names(weights) %in% regularized] <- 0
#' 
#' regsem <- elasticNet(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   weights = weights, 
#'   lambdas = seq(0,1,length.out = 20),
#'   alpha = seq(0,1,length.out = 5))
#' 
#' # elements of regsem can be accessed with the @ operator:
#' regsem@parameters[1,]
#' 
#' # use cross-validation to find the best parameters
#' cv <- cv4elasticNet(
#'   regularizedSEM = regsem, 
#'   k = 5
#' )
#' cv
#' 
#' #### Advanced ###
#' # Switching the optimizer # 
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' regsemGlmnet <- elasticNet(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   weights = weights, 
#'   lambdas = seq(0,1,length.out = 20),
#'   alpha = seq(0,1,length.out = 5),
#'   method = "glmnet",
#'   control = controlGlmnet())
#' 
#' # Note: The results are basically identical:
#' regsemGlmnet@parameters - regsem@parameters 
#' @export
elasticNet <- function(lavaanModel,
                       weights,
                       lambdas = NULL,
                       nLambdas = NULL,
                       alphas,
                       method = "ista", 
                       control = controlIsta()){
  
  inputArguments <- as.list(environment())
  
  if(! method %in% c("ista", "glmnet")) 
    stop("Currently ony methods = 'ista' and methods = 'glmnet' are supported")
  
  if(is.null(lambdas) && is.null(nLambdas)) 
    stop("Specify lambdas or nLambdas")
  if(!is.null(lambdas) && !is.null(nLambdas))
    stop("Specify either lambdas or nLambdas, not both")
  if(!is.null(nLambdas) && any(alphas != 1))
    stop("nLambdas only supported for lasso type penalty (alpha = 1).")
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet")
  
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
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
  }else if(any(startingValues == "start")){
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "start",
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
  }else if(is.numeric(startingValues)){
    
    if(!all(names(startingValues) %in% names(lessSEM::getLavaanParameters(lavaanModel)))) stop("Parameter names of startingValues do not match those of the lavaan object. See lessSEM::getLavaanParameters(lavaanModel).")
    SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                transformVariances = TRUE,
                                whichPars = "start", 
                                fit = FALSE,
                                addMeans = control$addMeans, 
                                activeSet = control$activeSet)
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
  
  #### glmnet requires an initial Hessian ####
  if(method == "glmnet"){
    initialHessian <- control$initialHessian
    if(is.matrix(initialHessian) && nrow(initialHessian) == length(rawParameters) && ncol(initialHessian) == length(rawParameters)){
      
      if(!all(rownames(initialHessian) %in% names(lessSEM::getLavaanParameters(lavaanModel))) ||
         !all(colnames(initialHessian) %in% names(lessSEM::getLavaanParameters(lavaanModel)))
      ) stop("initialHessian must have the parameter names as rownames and colnames. See lessSEM::getLavaanParameters(lavaanModel).")
      
    }else if(any(initialHessian == "compute")){
      
      initialHessian <- lessSEM:::getHessian(SEM = SEM, raw = TRUE)
      
    }else if(any(initialHessian == "scoreBased")){
      
      scores <- lessSEM:::getScores(SEM = SEM, raw = TRUE)
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
      stop("Invalid initialHessian passed to glmnet See ?controlGlmnet for more information.")
    }
    
    control$initialHessian <- initialHessian
  }
  
  #### prepare regularized model object ####
  if(method == "glmnet"){
    
    controlIntern <- list(
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
    
    regularizedModel <- new(glmnetEnet, 
                            weights, 
                            controlIntern)
    
  }else if(method == "ista"){
    
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
    
    regularizedModel <- new(istaEnet, 
                            weights, 
                            controlIntern)
  }
  
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
    
    inputArguments$lambdas = lambdas
    
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
  
  if(method == "glmnet" && control$saveHessian){
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
    
    if(method == "glmnet" && control$saveHessian) 
      Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning("Fit for lambda = ",
              lambda, "alpha = ", 
              alpha,
              " resulted in Error!")
      
      SEM <- lessSEM:::SEMFromLavaan(lavaanModel = lavaanModel, 
                                  transformVariances = TRUE,
                                  whichPars = startingValues,
                                  addMeans = control$addMeans)
      
      if(method == "glmnet"){
        regularizedModel$setHessian(controlIntern$initialHessian)
      }
      
    }else{
      
      if(method == "glmnet"){
        if(control$saveHessian) Hessians$Hessian[[it]] <- result$Hessian
        
        # set Hessian for next iteration
        regularizedModel$setHessian(result$Hessian)
        
      }
      
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