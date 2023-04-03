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
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param tau parameters below threshold tau will be seen as zeroed
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @returns Model of class regularizedSEM

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
#' # Regularization:
#' 
#' lsem <- smoothLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   epsilon = 1e-10,
#'   tau = 1e-4,
#'   lambdas = seq(0,1,length.out = 50))
#' 
#' # use the plot-function to plot the regularized parameters:
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' 
#' # AIC and BIC:
#' AIC(lsem)
#' BIC(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' coef(lsem, criterion = "BIC")
#' @export
smoothLasso <- function(lavaanModel,
                        regularized,
                        lambdas,
                        epsilon,
                        tau,
                        modifyModel = lessSEM::modifyModel(),
                        control = lessSEM::controlBFGS()){
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 1)
  
  result <- .regularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "lasso",
    weights = regularized,
    tuningParameters = tuningParameters,
    epsilon = epsilon, 
    tau = tau,
    modifyModel = modifyModel,
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
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param tau parameters below threshold tau will be seen as zeroed
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @returns Model of class regularizedSEM

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
#' # Regularization:
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
#' lsem <- smoothAdaptiveLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   regularized = regularized,
#'   weights = weights,
#'   epsilon = 1e-10,
#'   tau = 1e-4,
#'   lambdas = seq(0,1,length.out = 50))
#' 
#' # use the plot-function to plot the regularized parameters:
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' 
#' # AIC and BIC:
#' AIC(lsem)
#' BIC(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' coef(lsem, criterion = "BIC")
#' @export
smoothAdaptiveLasso <- function(lavaanModel,
                                regularized,
                                weights = NULL,
                                lambdas,
                                epsilon,
                                tau,
                                modifyModel = lessSEM::modifyModel(),
                                control = lessSEM::controlBFGS()){
  
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 1)
  
  if(is.null(weights)) weights <- regularized
  
  result <- .regularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "adaptiveLasso",
    weights = weights,
    tuningParameters = tuningParameters,
    epsilon = epsilon, 
    tau = tau,
    modifyModel = modifyModel,
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @returns Model of class regularizedSEM
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
#' # Regularization:
#' 
#' # names of the regularized parameters:
#' regularized = paste0("l", 6:15)
#' 
#' lsem <- ridgeBfgs(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   regularized = regularized,
#'   lambdas = seq(0,1,length.out = 50))
#' 
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' @export
ridgeBfgs <- function(lavaanModel,
                      regularized,
                      lambdas = NULL,
                      modifyModel = lessSEM::modifyModel(),
                      control = lessSEM::controlBFGS()){
  
  result <- .regularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "ridge",
    weights = regularized,
    tuningParameters = data.frame(lambda = lambdas,
                                  alpha = 0),
    epsilon = 0, # ridge is already smooth
    tau = 0,
    modifyModel = modifyModel,
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
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' between 0 and 1. 0 = ridge, 1 = lasso.
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param tau parameters below threshold tau will be seen as zeroed
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @returns Model of class regularizedSEM

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
#' # Regularization:
#' 
#' # names of the regularized parameters:
#' regularized = paste0("l", 6:15)
#' 
#' lsem <- smoothElasticNet(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   regularized = regularized,
#'   epsilon = 1e-10,
#'   tau = 1e-4,
#'   lambdas = seq(0,1,length.out = 5),
#'   alphas = seq(0,1,length.out = 3))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' @export
smoothElasticNet <- function(lavaanModel,
                             regularized,
                             lambdas = NULL,
                             nLambdas = NULL,
                             alphas,
                             epsilon, 
                             tau,
                             modifyModel = lessSEM::modifyModel(),
                             control = lessSEM::controlBFGS()){
  
  if(any(alphas < 0) || any(alphas > 1)) 
    stop("alpha must be between 0 and 1.")
  
  result <- .regularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "elasticNet",
    weights = regularized,
    tuningParameters = expand.grid(lambda = lambdas,
                                   alpha = alphas),
    epsilon = epsilon, 
    tau = tau,
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}