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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
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
                  modifyModel = lessSEM::modifyModel(),
                  control = controlIsta()){
  
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas)
  }else{
    tuningParameters <- data.frame(lambdas = lambdas,
                                   alpha = 1)
  }
  
  result <- regularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "lasso",
    weights = regularized,
    tuningParameters = tuningParameters,
    method = method, 
    modifyModel = modifyModel,
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
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
                          modifyModel = lessSEM::modifyModel(),
                          control = lessSEM::controlIsta()){
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas)
  }else{
    tuningParameters <- data.frame(lambdas = lambdas,
                                   alpha = 1)
  }
  
  if(is.null(weights)) weights <- regularized
  
  result <- regularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "adaptiveLasso",
    weights = weights,
    tuningParameters = tuningParameters,
    method = method, 
    modifyModel = modifyModel,
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
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
                  lambdas,
                  method = "ista", 
                  modifyModel = lessSEM::modifyModel(),
                  control = controlIsta()){
  
  
  result <- regularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "ridge",
    weights = regularized,
    tuningParameters = data.frame(lambda = lambdas,
                                   alpha = 0),
    method = method, 
    modifyModel = modifyModel,
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
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
                       regularized,
                       lambdas,
                       alphas,
                       method = "ista", 
                       modifyModel = lessSEM::modifyModel(),
                       control = controlIsta()){
  
  if(any(alphas < 0) || any(alphas > 1)) 
    stop("alpha must be between 0 and 1.")
  
  result <- regularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "elasticNet",
    weights = regularized,
    tuningParameters = expand.grid(lambda = lambdas,
                                   alpha = alphas),
    method = method, 
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}