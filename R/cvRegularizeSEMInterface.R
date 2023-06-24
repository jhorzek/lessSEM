#' cvLasso
#' 
#' Implements cross-validated lasso regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda |x_j|}
#' Lasso regularization will set parameters to zero if \eqn{\lambda} is large enough
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
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
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
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
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
#' lsem <- cvLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1),
#'   k = 5, # number of cross-validation folds
#'   standardize = TRUE) # automatic standardization
#' 
#' # use the plot-function to plot the cross-validation fit:
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # The best parameters can also be extracted with:
#' estimates(lsem)
#' @export
cvLasso <- function(lavaanModel,
                    regularized,
                    lambdas,
                    k = 5,
                    standardize = FALSE,
                    returnSubsetParameters = FALSE,
                    method = "glmnet", 
                    modifyModel = lessSEM::modifyModel(),
                    control = lessSEM::controlGlmnet()){
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 1)
  
  result <- .cvRegularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "lasso",
    k = k,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    weights = regularized,
    tuningParameters = tuningParameters,
    method = method, 
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}

#' cvAdaptiveLasso
#' 
#' Implements cross-validated adaptive lasso regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = p( x_j) = \frac{1}{w_j}\lambda| x_j|}
#' Adaptive lasso regularization will set parameters to zero if \eqn{\lambda} 
#' is large enough.
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currenlty,
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
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
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
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
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
#' lsem <- cvAdaptiveLasso(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,.1))
#' 
#' # use the plot-function to plot the cross-validation fit
#' plot(lsem)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # The best parameters can also be extracted with:
#' estimates(lsem)
#' @export
cvAdaptiveLasso <- function(lavaanModel,
                            regularized,
                            weights = NULL,
                            lambdas,
                            k = 5,
                            standardize = FALSE,
                            returnSubsetParameters = FALSE,
                            method = "glmnet", 
                            modifyModel = lessSEM::modifyModel(),
                            control = lessSEM::controlGlmnet()){
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 1)
  
  if(is.null(weights)) weights <- regularized
  
  result <- .cvRegularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "adaptiveLasso",
    weights = weights,
    k = k,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    tuningParameters = tuningParameters,
    method = method, 
    modifyModel = modifyModel,
    control = control
  )
  
  return(result)
  
}

#' cvRidge
#' 
#' Implements ridge regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda x_j^2}
#' Note that ridge regularization will not set any of the parameters to zero
#' but result in a shrinkage towards zero. 
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currenlty,
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
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
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
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
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
#' lsem <- cvRidge(
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
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' @export
cvRidge <- function(lavaanModel,
                    regularized,
                    lambdas,
                    k = 5,
                    standardize = FALSE,
                    returnSubsetParameters = FALSE,
                    method = "glmnet", 
                    modifyModel = lessSEM::modifyModel(),
                    control = lessSEM::controlGlmnet()){
  
  
  result <- .cvRegularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "ridge",
    weights = regularized,
    k = k,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    tuningParameters = data.frame(lambda = lambdas,
                                  alpha = 0),
    method = method, 
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}

#' cvElasticNet
#' 
#' Implements elastic net regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \alpha\lambda| x_j| + (1-\alpha)\lambda x_j^2}
#' Note that the elastic net combines ridge and lasso regularization. If \eqn{\alpha = 0}, 
#' the elastic net reduces to ridge regularization. If \eqn{\alpha = 1} it reduces
#' to lasso regularization. In between, elastic net is a compromise between the shrinkage of
#' the lasso and the ridge penalty. 
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currenlty,
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
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#' 
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' between 0 and 1. 0 = ridge, 1 = lasso.
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures.
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
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
#' lsem <- cvElasticNet(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 5),
#'   alphas = seq(0,1,length.out = 3))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # optional: plotting the cross-validation fit requires installation of plotly
#' # plot(lsem)
#' @export
cvElasticNet <- function(lavaanModel,
                         regularized,
                         lambdas,
                         alphas,
                         k = 5,
                         standardize = FALSE,
                         returnSubsetParameters = FALSE,
                         method = "glmnet", 
                         modifyModel = lessSEM::modifyModel(),
                         control = lessSEM::controlGlmnet()){
  
  if(any(alphas < 0) || any(alphas > 1)) 
    stop("alpha must be between 0 and 1.")
  
  result <- .cvRegularizeSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "elasticNet",
    weights = regularized,
    k = k,
    standardize = standardize,
    returnSubsetParameters = returnSubsetParameters,
    tuningParameters = expand.grid(lambda = lambdas,
                                   alpha = alphas),
    method = method, 
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}


#' cvCappedL1
#' 
#' Implements cappedL1 regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \min(| x_j|, \theta)}
#' where \eqn{\theta > 0}. The cappedL1 penalty is identical to the lasso for 
#' parameters which are below \eqn{\theta} and identical to a constant for parameters
#' above \eqn{\theta}. As adding a constant to the fitting function will not change its
#' minimum, larger parameters can stay unregularized while smaller ones are set to zero.
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currenlty,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' CappedL1 regularization:
#' 
#' * Zhang, T. (2010). Analysis of Multi-stage Convex Relaxation for Sparse Regularization. 
#' Journal of Machine Learning Research, 11, 1081–1107.
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta function. See ?controlIsta
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
#' lsem <- cvCappedL1(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 5),
#'   thetas = seq(0.01,2,length.out = 3))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # optional: plotting the cross-validation fit requires installation of plotly
#' # plot(lsem)
#' @export
cvCappedL1 <- function(lavaanModel,
                       regularized,
                       lambdas,
                       thetas,
                       k = 5,
                       standardize = FALSE,
                       returnSubsetParameters = FALSE,
                       modifyModel = lessSEM::modifyModel(),
                       method = "glmnet",
                       control = lessSEM::controlGlmnet()){
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .cvRegularizeSEMInternal(lavaanModel = lavaanModel, 
                                     penalty = "cappedL1", 
                                     weights = regularized,
                                     k = k,
                                     standardize = standardize,
                                     returnSubsetParameters = returnSubsetParameters,
                                     tuningParameters = expand.grid(lambda = lambdas, 
                                                                    theta = thetas,
                                                                    alpha = 1), 
                                     method = method, 
                                     modifyModel = modifyModel, 
                                     control = control
  )
  
  return(result)
  
}

#' cvLsp
#' 
#' Implements lsp regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \log(1 + |x_j|/\theta)}
#' where \eqn{\theta > 0}. 
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currenlty,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' lsp regularization:
#' 
#' * Candès, E. J., Wakin, M. B., & Boyd, S. P. (2008). Enhancing Sparsity by 
#' Reweighted l1 Minimization. Journal of Fourier Analysis and Applications, 14(5–6), 
#' 877–905. https://doi.org/10.1007/s00041-008-9045-x
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta function. See ?controlIsta
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
#' lsem <- cvLsp(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 5),
#'   thetas = seq(0.01,2,length.out = 3))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # optional: plotting the cross-validation fit requires installation of plotly
#' # plot(lsem)
#' @export
cvLsp <- function(lavaanModel,
                  regularized,
                  lambdas,
                  thetas,
                  k = 5,
                  standardize = FALSE,
                  returnSubsetParameters = FALSE,
                  modifyModel = lessSEM::modifyModel(),
                  method = "glmnet",
                  control = lessSEM::controlGlmnet()){
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .cvRegularizeSEMInternal(lavaanModel = lavaanModel, 
                                     penalty = "lsp", 
                                     weights = regularized,
                                     k = k,
                                     standardize = standardize,
                                     returnSubsetParameters = returnSubsetParameters,
                                     tuningParameters = expand.grid(lambda = lambdas, 
                                                                    theta = thetas), 
                                     method = method, 
                                     modifyModel = modifyModel, 
                                     control = control
  )
  
  return(result)
  
}

#' cvMcp
#' 
#' Implements mcp regularization for structural equation models.
#' The penalty function is given by:
#' \ifelse{html}{\deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| - x_j^2/(2\theta) & \text{if } |x_j| \leq \theta\lambda\\
#' \theta\lambda^2/2 & \text{if } |x_j| > \lambda\theta
#' \end{cases}} where \eqn{\theta > 0}.}{
#' Equation Omitted in Pdf Documentation.}
#' 
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currenlty,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' mcp regularization:
#' 
#' * Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax concave penalty. 
#' The Annals of Statistics, 38(2), 894–942. https://doi.org/10.1214/09-AOS729
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta function. See ?controlIsta
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
#' lsem <- cvMcp(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 5),
#'   thetas = seq(0.01,2,length.out = 3))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # optional: plotting the cross-validation fit requires installation of plotly
#' # plot(lsem)
#' @export
cvMcp <- function(lavaanModel,
                  regularized,
                  lambdas,
                  thetas,
                  k = 5,
                  standardize = FALSE,
                  returnSubsetParameters = FALSE,
                  modifyModel = lessSEM::modifyModel(),
                  method = "ista",
                  control = lessSEM::controlIsta()){
  
  if(any(thetas <= 0)) 
    stop("Theta must be > 0")
  if(any(thetas <= 1) & method == "glmnet") 
    warning("thetas is typically > 1. Note that glmnet may run into issues with small theta.")
  result <- .cvRegularizeSEMInternal(lavaanModel = lavaanModel, 
                                     penalty = "mcp", 
                                     weights = regularized,
                                     k = k,
                                     standardize = standardize,
                                     returnSubsetParameters = returnSubsetParameters,
                                     tuningParameters = expand.grid(lambda = lambdas, 
                                                                    theta = thetas), 
                                     method = method, 
                                     modifyModel = modifyModel, 
                                     control = control
  )
  
  return(result)
  
}

#' cvScad
#' 
#' Implements scad regularization for structural equation models.
#' The penalty function is given by:
#' \ifelse{html}{
#' \deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| & \text{if } |x_j| \leq \theta\\
#' \frac{-x_j^2 + 2\theta\lambda |x_j| - \lambda^2}{2(\theta -1)} & 
#' \text{if } \lambda < |x_j| \leq \lambda\theta \\
#' (\theta + 1) \lambda^2/2 & \text{if } |x_j| \geq \theta\lambda\\
#' \end{cases}}
#' where \eqn{\theta > 2}.}{
#' Equation Omitted in Pdf Documentation.
#' }
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currenlty,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' scad regularization:
#' 
#' * Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized 
#' likelihood and its oracle properties. Journal of the American Statistical Association, 
#' 96(456), 1348–1360. https://doi.org/10.1198/016214501753382273
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#' 
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
#' * Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param lavaanModel model of class lavaan 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param k the number of cross-validation folds. Alternatively, you can pass
#' a matrix with booleans (TRUE, FALSE) which indicates for each person which subset
#' it belongs to. See ?lessSEM::createSubsets for an example of how this matrix should look like.
#' @param standardize Standardizing your data prior to the analysis can undermine the cross-
#' validation. Set standardize=TRUE to automatically standardize the data.
#' @param returnSubsetParameters set to TRUE to return the parameters for each training set
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista and glmnet. 
#' With ista, the control argument can be used to switch to related procedures.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta function. See ?controlIsta
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
#' lsem <- cvScad(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 3),
#'   thetas = seq(2.01,5,length.out = 3))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' 
#' # optional: plotting the cross-validation fit requires installation of plotly
#' # plot(lsem)
#' @export
cvScad <- function(lavaanModel,
                   regularized,
                   lambdas,
                   thetas,
                   k = 5,
                   standardize = FALSE,
                   returnSubsetParameters = FALSE,
                   modifyModel = lessSEM::modifyModel(),
                   method = "glmnet",
                   control = lessSEM::controlGlmnet()){
  
  if(any(thetas <= 2)) stop("Theta must be > 2")
  
  result <- .cvRegularizeSEMInternal(lavaanModel = lavaanModel, 
                                     penalty = "scad", 
                                     weights = regularized,
                                     k = k,
                                     standardize = standardize,
                                     returnSubsetParameters = returnSubsetParameters,
                                     tuningParameters = expand.grid(lambda = lambdas, 
                                                                    theta = thetas), 
                                     method = method, 
                                     modifyModel = modifyModel, 
                                     control = control
  )
  
  return(result)
}