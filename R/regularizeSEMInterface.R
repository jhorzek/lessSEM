#' lasso
#' 
#' Implements lasso regularization for structural equation models.
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
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param reverse if set to TRUE and nLambdas is used, lessSEM will start with the
#' largest lambda and gradually decrease lambda. Otherwise, lessSEM will start with
#' the smallest lambda and gradually increase it.
#' @param curve Allows for unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). 
#' If curve is close to 1 all lambda values will be equally spaced, if curve is large 
#' lambda values will be more concentrated close to 0. See ?lessSEM::curveLambda for more information.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
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
#' lsem <- lasso(
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
#' lsem@parameters[1,]
#' 
#' # fit Measures:
#' fitIndices(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' # or
#' estimates(lsem, criterion = "AIC") 
#' 
#' #### Advanced ###
#' # Switching the optimizer # 
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' lsemIsta <- lasso(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   nLambdas = 50,
#'   method = "ista",
#'   control = controlIsta())
#' 
#' # Note: The results are basically identical:
#' lsemIsta@parameters - lsem@parameters
#' @export
lasso <- function(lavaanModel,
                  regularized,
                  lambdas = NULL,
                  nLambdas = NULL,
                  reverse = TRUE,
                  curve = 1,
                  method = "glmnet", 
                  modifyModel = lessSEM::modifyModel(),
                  control = lessSEM::controlGlmnet()){
  
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas,
                                   reverse = reverse,
                                   curve = curve)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1)
  }
  
  result <- .regularizeSEMInternal(
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
#' Implements adaptive lasso regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = p( x_j) = \frac{1}{w_j}\lambda| x_j|}
#' Adaptive lasso regularization will set parameters to zero if \eqn{\lambda} 
#' is large enough.
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
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param reverse if set to TRUE and nLambdas is used, lessSEM will start with the
#' largest lambda and gradually decrease lambda. Otherwise, lessSEM will start with
#' the smallest lambda and gradually increase it.
#' @param curve Allows for unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). 
#' If curve is close to 1 all lambda values will be equally spaced, if curve is large 
#' lambda values will be more concentrated close to 0. See ?lessSEM::curveLambda for more information.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
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
#' lsem <- adaptiveLasso(
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
#' lsem@parameters[1,]
#' 
#' # fit Measures:
#' fitIndices(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' # or
#' estimates(lsem, criterion = "AIC")
#' 
#' #### Advanced ###
#' # Switching the optimizer #
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' lsemIsta <- adaptiveLasso(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   nLambdas = 50,
#'   method = "ista",
#'   control = controlIsta())
#' 
#' # Note: The results are basically identical:
#' lsemIsta@parameters - lsem@parameters
#' @export
adaptiveLasso <- function(lavaanModel,
                          regularized,
                          weights = NULL,
                          lambdas = NULL,
                          nLambdas = NULL,
                          reverse = TRUE,
                          curve = 1,
                          method = "glmnet", 
                          modifyModel = lessSEM::modifyModel(),
                          control = lessSEM::controlGlmnet()){
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas,
                                   reverse = reverse,
                                   curve = curve)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1)
  }
  
  if(is.null(weights)) weights <- regularized
  
  result <- .regularizeSEMInternal(
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
#' Implements ridge regularization for structural equation models.
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
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
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
#' lsem <- ridge(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20))
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
#' 
#' #### Advanced ###
#' # Switching the optimizer #
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' lsemIsta <- ridge(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20),
#'   method = "ista",
#'   control = controlIsta())
#' 
#' # Note: The results are basically identical:
#' lsemIsta@parameters - lsem@parameters
#' @export
ridge <- function(lavaanModel,
                  regularized,
                  lambdas,
                  method = "glmnet", 
                  modifyModel = lessSEM::modifyModel(),
                  control = lessSEM::controlGlmnet()){
  
  
  result <- .regularizeSEMInternal(
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
#' Implements elastic net regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \alpha\lambda| x_j| + (1-\alpha)\lambda x_j^2}
#' Note that the elastic net combines ridge and lasso regularization. If \eqn{\alpha = 0}, 
#' the elastic net reduces to ridge regularization. If \eqn{\alpha = 1} it reduces
#' to lasso regularization. In between, elastic net is a compromise between the shrinkage of
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
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. This element is generated with 
#' the lessSEM::controlIsta() and controlGlmnet() functions.
#' @returns Model of class regularizedSEM
#'
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
#' lsem <- elasticNet(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 5),
#'   alphas = seq(0,1,length.out = 3))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' 
#' # optional: plotting the paths requires installation of plotly
#' # plot(lsem)
#' 
#' #### Advanced ###
#' # Switching the optimizer #
#' # Use the "method" argument to switch the optimizer. The control argument
#' # must also be changed to the corresponding function:
#' lsemIsta <- elasticNet(
#'   lavaanModel = lavaanModel,
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 5),
#'   alphas = seq(0,1,length.out = 3),
#'   method = "ista",
#'   control = controlIsta())
#' 
#' # Note: The results are basically identical:
#' lsemIsta@parameters - lsem@parameters
#' @export
elasticNet <- function(lavaanModel,
                       regularized,
                       lambdas,
                       alphas,
                       method = "glmnet", 
                       modifyModel = lessSEM::modifyModel(),
                       control = lessSEM::controlGlmnet()){
  
  if(any(alphas < 0) || any(alphas > 1)) 
    stop("alpha must be between 0 and 1.")
  
  result <- .regularizeSEMInternal(
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


#' cappedL1
#' 
#' Implements cappedL1 regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \min(| x_j|, \theta)}
#' where \eqn{\theta > 0}. The cappedL1 penalty is identical to the lasso for 
#' parameters which are below \eqn{\theta} and identical to a constant for parameters
#' above \eqn{\theta}. As adding a constant to the fitting function will not change its
#' minimum, larger parameters can stay unregularized while smaller ones are set to zero.
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta (see ?controlIsta)
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
#' lsem <- cappedL1(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20),
#'   thetas = seq(0.01,2,length.out = 5))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' 
#' # fit Measures:
#' fitIndices(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' # or
#' estimates(lsem, criterion = "AIC")
#' 
#' # optional: plotting the paths requires installation of plotly
#' # plot(lsem)
#' @export
cappedL1 <- function(lavaanModel,
                     regularized,
                     lambdas,
                     thetas,
                     modifyModel = lessSEM::modifyModel(),
                     method = "glmnet",
                     control = lessSEM::controlGlmnet()){
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .regularizeSEMInternal(lavaanModel = lavaanModel, 
                                           penalty = "cappedL1", 
                                           weights = regularized,
                                           tuningParameters = expand.grid(lambda = lambdas, 
                                                                          theta = thetas,
                                                                          alpha = 1), 
                                           method = method, 
                                           modifyModel = modifyModel, 
                                           control = control
  )
  
  return(result)
  
}

#' lsp
#' 
#' Implements lsp regularization for structural equation models.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \log(1 + |x_j|/\theta)}
#' where \eqn{\theta > 0}. 
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta (see ?controlIsta)
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
#' lsem <- lsp(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20),
#'   thetas = seq(0.01,2,length.out = 5))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' 
#' # fit Measures:
#' fitIndices(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' # or
#' estimates(lsem, criterion = "AIC")
#' 
#' # optional: plotting the paths requires installation of plotly
#' # plot(lsem)
#' @export
lsp <- function(lavaanModel,
                regularized,
                lambdas,
                thetas,
                modifyModel = lessSEM::modifyModel(),
                method = "glmnet",
                control = lessSEM::controlGlmnet()){
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .regularizeSEMInternal(lavaanModel = lavaanModel, 
                                           penalty = "lsp", 
                                           weights = regularized,
                                           tuningParameters = expand.grid(lambda = lambdas, 
                                                                          theta = thetas), 
                                           method = method, 
                                           modifyModel = modifyModel, 
                                           control = control
  )
  
  return(result)
  
}

#' mcp
#' 
#' Implements mcp regularization for structural equation models.
#' The penalty function is given by:
#' \ifelse{html}{\deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| - x_j^2/(2\theta) & \text{if } |x_j| \leq \theta\lambda\\
#' \theta\lambda^2/2 & \text{if } |x_j| > \lambda\theta
#' \end{cases}} where \eqn{\theta > 1}.}{
#' Equation Omitted in Pdf Documentation.}
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well.
#' 
#' In our experience, the glmnet optimizer can run in issues with the mcp penalty.
#' Therefor, we default to using ista.
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta (see ?controlIsta)
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
#' lsem <- mcp(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20),
#'   thetas = seq(0.01,2,length.out = 5))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' 
#' # fit Measures:
#' fitIndices(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' # or
#' estimates(lsem, criterion = "AIC")
#' 
#' # optional: plotting the paths requires installation of plotly
#' # plot(lsem)
#' @export
mcp <- function(lavaanModel,
                regularized,
                lambdas,
                thetas,
                modifyModel = lessSEM::modifyModel(),
                method = "ista",
                control = lessSEM::controlIsta()){
  
  if(any(thetas <= 0)) 
    stop("Theta must be > 0")
  if(any(thetas <= 1) & method == "glmnet") 
    warning("thetas is typically > 1. Note that glmnet may run into issues with small theta.")
  result <- .regularizeSEMInternal(lavaanModel = lavaanModel, 
                                           penalty = "mcp", 
                                           weights = regularized,
                                           tuningParameters = expand.grid(lambda = lambdas, 
                                                                          theta = thetas), 
                                           method = method, 
                                           modifyModel = modifyModel, 
                                           control = control
  )
  
  return(result)
  
}

#' scad
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
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
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
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta (see ?controlIsta)
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
#' lsem <- scad(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel,
#'   # names of the regularized parameters:
#'   regularized = paste0("l", 6:15),
#'   lambdas = seq(0,1,length.out = 20),
#'   thetas = seq(2.01,5,length.out = 5))
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # if you are only interested in the estimates and not the tuning parameters, use
#' coef(lsem)@estimates
#' # or
#' estimates(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters[1,]
#' 
#' # fit Measures:
#' fitIndices(lsem)
#' 
#' # The best parameters can also be extracted with:
#' coef(lsem, criterion = "AIC")
#' # or
#' estimates(lsem, criterion = "AIC")
#' 
#' # optional: plotting the paths requires installation of plotly
#' # plot(lsem)
#' @export
scad <- function(lavaanModel,
                 regularized,
                 lambdas,
                 thetas,
                 modifyModel = lessSEM::modifyModel(),
                 method = "glmnet",
                 control = lessSEM::controlGlmnet()){
  
  if(any(thetas <= 2)) stop("Theta must be > 2")
  
  result <- .regularizeSEMInternal(lavaanModel = lavaanModel, 
                                           penalty = "scad", 
                                           weights = regularized,
                                           tuningParameters = expand.grid(lambda = lambdas, 
                                                                          theta = thetas), 
                                           method = method, 
                                           modifyModel = modifyModel, 
                                           control = control
  )
  
  return(result)
}