#' cvSmoothLasso
#' 
#' Implements cross-validated smooth lasso regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \sqrt{(x_j + \epsilon)^2}}
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' Lasso regularization:
#' 
#' * Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical 
#' Society. Series B (Methodological), 58(1), 267–288.
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @returns model of class cvRegularizedSEM

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
#' lsem <- cvSmoothLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1),
#'   k = 5, # number of cross-validation folds
#'   epsilon = 1e-8,
#'   standardize = TRUE) # automatic standardization
#' 
#' # use the plot-function to plot the cross-validation fit:
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem)
#' @export
cvSmoothLasso <- function(lavaanModel,
                          regularized,
                          lambdas,
                          epsilon,
                          k = 5,
                          standardize = FALSE,
                          returnSubsetParameters = FALSE,
                          modifyModel = lessSEM::modifyModel(),
                          control = lessSEM::controlBFGS()){
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 1)
  
  result <- .cvRegularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "lasso",
    k = k,
    epsilon = epsilon,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    weights = regularized,
    tuningParameters = tuningParameters,
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}

#' cvSmoothAdaptiveLasso
#' 
#' Implements cross-validated smooth adaptive lasso regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = p( x_j) = \frac{1}{w_j}\lambda\sqrt{(x_j + \epsilon)^2}}
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' Adaptive lasso regularization:
#' 
#' * Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American Statistical Association, 
#' 101(476), 1418–1429. https://doi.org/10.1198/016214506000000735
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
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
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS for more details.
#' @returns model of class cvRegularizedSEM 

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
#' lsem <- cvSmoothAdaptiveLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1),
#'   epsilon = 1e-8)
#' 
#' # use the plot-function to plot the cross-validation fit
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem)
#' @export
cvSmoothAdaptiveLasso <- function(lavaanModel,
                            regularized,
                            weights = NULL,
                            lambdas,
                            epsilon,
                            k = 5,
                            standardize = FALSE,
                            returnSubsetParameters = FALSE,
                            modifyModel = lessSEM::modifyModel(),
                            control = lessSEM::controlBFGS()){
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 1)
  
  if(is.null(weights)) weights <- regularized
  
  result <- .cvRegularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "adaptiveLasso",
    weights = weights,
    epsilon = epsilon,
    k = k,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    tuningParameters = tuningParameters,
    modifyModel = modifyModel,
    control = control
  )
  
  return(result)
  
}

#' cvRidgeBfgs
#' 
#' Implements cross-validated ridge regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda x_j^2}
#' Note that ridge regularization will not set any of the parameters to zero
#' but result in a shrinkage towards zero. 
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' Ridge regularization:
#' 
#' * Hoerl, A. E., & Kennard, R. W. (1970). Ridge Regression: Biased Estimation 
#' for Nonorthogonal Problems. Technometrics, 12(1), 55–67. 
#' https://doi.org/10.1080/00401706.1970.10488634
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS
#' for more details.
#' @returns model of class cvRegularizedSEM

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
#' lsem <- cvRidgeBfgs(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20))
#' 
#' # use the plot-function to plot the cross-validation fit:
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' @export
cvRidgeBfgs <- function(lavaanModel,
                    regularized,
                    lambdas,
                    k = 5,
                    standardize = FALSE,
                    returnSubsetParameters = FALSE,
                    modifyModel = lessSEM::modifyModel(),
                    control = lessSEM::controlBFGS()){
  
  
  result <- .cvRegularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "ridge",
    weights = regularized,
    epsilon = 0,
    k = k,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    tuningParameters = data.frame(lambda = lambdas,
                                  alpha = 0),
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}

#' cvSmoothElasticNet
#' 
#' Implements cross-validated  smooth elastic net regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \alpha\lambda\sqrt{(x_j + \epsilon)^2} + (1-\alpha)\lambda x_j^2}
#' Note that the smooth elastic net combines ridge and smooth lasso regularization. If \eqn{\alpha = 0}, 
#' the elastic net reduces to ridge regularization. If \eqn{\alpha = 1} it reduces
#' to smooth lasso regularization. In between, elastic net is a compromise between the shrinkage of
#' the lasso and the ridge penalty. 
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' Elastic net regularization:
#' 
#' * Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. 
#' Journal of the Royal Statistical Society: Series B, 67(2), 301–320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' between 0 and 1. 0 = ridge, 1 = lasso.
#' @param epsilon epsilon > 0; controls the smoothness of the approximation. Larger values = smoother 
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlBFGS function. See ?controlBFGS
#' for more details.
#' @returns model of class cvRegularizedSEM
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
#' lsem <- cvSmoothElasticNet(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   epsilon = 1e-8,
#'   lambdas = seq(0,1,length.out = 5),
#'   alphas = .3)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # optional: plotting the cross-validation fit requires installation of plotly
#' # plot(lsem)
#' @export
cvSmoothElasticNet <- function(lavaanModel,
                         regularized,
                         lambdas,
                         alphas,
                         epsilon,
                         k = 5,
                         standardize = FALSE,
                         returnSubsetParameters = FALSE,
                         modifyModel = lessSEM::modifyModel(),
                         control = lessSEM::controlBFGS()){
  
  if(any(alphas < 0) || any(alphas > 1)) 
    stop("alpha must be between 0 and 1.")
  
  result <- .cvRegularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "elasticNet",
    weights = regularized,
    epsilon = epsilon,
    k = k,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    tuningParameters = expand.grid(lambda = lambdas,
                                   alpha = alphas),
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}
