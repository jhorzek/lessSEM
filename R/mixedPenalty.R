#' mixedPenalty
#' 
#' Provides possibility to impose different penalties on different parameters.
#' 
#' The `mixedPenalty` function allows you to add multiple penalties to a single model.
#' For instance, you may want to regularize both loadings and regressions in a SEM.
#' In this case, using the same penalty (e.g., lasso) for both types of penalties may
#' actually not be what you want to use because the penalty function is sensitive to
#' the scales of the parameters. Instead, you may want to use two separate lasso
#' penalties for loadings and regressions. Similarly, separate penalties for 
#' different parameters have, for instance, been proposed in multi-group models
#' (Geminiani et al., 2021).
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
#' most standard SEM are supported. \pkg{lessSEM} also provides full information
#' maximum likelihood for missing data. To use this functionality,
#' fit your \pkg{lavaan} model with the argument `sem(..., missing = 'ml')`. 
#' \pkg{lessSEM} will then automatically switch to full information maximum likelihood
#' as well. Models are fitted with the glmnet or ista optimizer. Note that the 
#' optimizers differ in which penalties they support. The following table provides
#' an overview:
#' 
#' | Penalty | Function | glmnet | ista |
#' | --- | ---- | ---- | ---- |
#' | lasso | addLasso | x | x |
#' | elastic net | addElasticNet | x* | - |
#' | cappedL1 | addCappedL1 | x | x |
#' | lsp | addLsp | x | x |
#' | scad | addScad | x | x |
#' | mcp | addMcp | x | x |
#' 
#' By default, glmnet will be used. Note that the elastic net penalty can
#' only be combined with other elastic net penalties.
#' 
#' Check vignette(topic = "Mixed-Penalties", package = "lessSEM") for more details.
#' 
#' Regularized SEM
#' 
#' * Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
#' * Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural 
#' Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
#'  
#' For more details on ISTA, see:
#' 
#' * Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' * Geminiani, E., Marra, G., & Moustaki, I. (2021). Single- and multiple-group penalized factor analysis: 
#' A trust-region algorithm approach with integrated automatic multiple tuning parameter selection. 
#' Psychometrika, 86(1), 65–95. https://doi.org/10.1007/s11336-021-09751-8
#' * Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#' * Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and 
#' Trends in Optimization, 1(3), 123–231.
#'  
#' @param lavaanModel model of class lavaan 
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param method which optimizer should be used? Currently supported are "glmnet" and "ista".
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
#' # In this example, we want to regularize the loadings l6-l10 
#' # independently of the loadings l11-15. This could, for instance,
#' # reflect that the items y6-y10 and y11-y15 may belong to different
#' # subscales. 
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add lasso penalty on loadings l6 - l10:
#'   addLasso(regularized = paste0("l", 6:10), 
#'            lambdas = seq(0,1,length.out = 4)) |>
#'   # add scad penalty on loadings l11 - l15:
#'   addScad(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,length.out = 3),
#'           thetas = 3.1) |>
#'   # fit the model:
#'   fit()
#' 
#' # elements of regularized can be accessed with the @ operator:
#' regularized@parameters[1,]
#' 
#' # AIC and BIC:
#' AIC(regularized)
#' BIC(regularized)
#' 
#' # The best parameters can also be extracted with:
#' coef(regularized, criterion = "AIC")
#' coef(regularized, criterion = "BIC")
#' 
#' # The tuningParameterConfiguration corresponds to the rows
#' # in the lambda, theta, and alpha matrices in regularized@tuningParamterConfigurations.
#' # Configuration 3, for example, is given by
#' regularized@tuningParameterConfigurations$lambda[3,]
#' regularized@tuningParameterConfigurations$theta[3,]
#' regularized@tuningParameterConfigurations$alpha[3,] 
#' # Note that lambda, theta, and alpha may correspond to tuning parameters
#' # of different penalties for different parameters (e.g., lambda for l6 is the lambda
#' # of the lasso penalty, while lambda for l12 is the lambda of the scad penalty).
#' 
#' @md
#' @export
mixedPenalty <- function(lavaanModel,
                         modifyModel = lessSEM::modifyModel(),
                         method = "glmnet", 
                         control = lessSEM::controlGlmnet()){
  notes <- c("Notes:")
  
  notes <- c(
    notes,
    paste0("Mixed penalties is a very new feature. Please note that there may still ",
                  "be bugs in the procedure. Use carefully!"))
  
  mixedPenalty <- list(
    lavaanModel = lavaanModel,
    method = method,
    modifyModel = modifyModel,
    control = control,
    penalties = list()
  )
  
  class(mixedPenalty) <- "mixedPenalty"
  
  return(mixedPenalty)
}


#' addCappedL1
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
#' @param mixedPenalty model of class mixedPenalty created with the mixedPenalty function (see ?mixedPenalty) 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @returns Model of class mixedPenalty. Use the fit() - function to fit the model
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
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addCappedL1(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,.1),
#'           thetas = 2.3) |>
#'   # fit the model:
#'   fit()
#' @export
addCappedL1 <- function(mixedPenalty,
                        regularized,
                        lambdas,
                        thetas){
  
  if(!is(mixedPenalty, "mixedPenalty"))
    stop("mixedPenalty must be of class mixedPenalty. ",
         "These models can be created with the regularize() function.")
  
  tps <- expand.grid(lambda = lambdas,
                     theta = thetas,
                     alpha = 1)
  
  mixedPenalty$penalties <- c(
    mixedPenalty$penalties,
    list(list(regularized = regularized, 
              tps = tps,
              penaltyType = .penaltyTypes("cappedL1"))
    )
  )
  
  return(mixedPenalty)
}

#' addLasso
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
#' @param mixedPenalty model of class mixedPenalty created with the mixedPenalty function (see ?mixedPenalty) 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param weights can be used to give different weights to the different parameters
#' @returns Model of class mixedPenalty. Use the fit() - function to fit the model
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
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addLasso(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,.1)) |>
#'   # fit the model:
#'   fit()
#' @export
addLasso <- function(mixedPenalty,
                     regularized,
                     weights = 1,
                     lambdas){
  
  if(!is(mixedPenalty, "mixedPenalty"))
    stop("mixedPenalty must be of class mixedPenalty. ",
         "These models can be created with the regularize() function.")
  
  if((length(weights) != 1) && (length(weights) != length(regularized)))
    stop("Weights argument must be of length 1 or of length(regularized).")
  
  tps <- expand.grid(lambda = lambdas,
                     theta = 0,
                     alpha = 1)
  
  mixedPenalty$penalties <- c(
    mixedPenalty$penalties,
    list(list(regularized = regularized,
              weight = weights,
              tps = tps,
              penaltyType = .penaltyTypes("lasso")))
  )
  
  return(mixedPenalty)
}

#' addLsp
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
#' @param mixedPenalty model of class mixedPenalty created with the mixedPenalty function (see ?mixedPenalty) 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @returns Model of class mixedPenalty. Use the fit() - function to fit the model
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
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addLsp(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,.1),
#'           thetas = 2.3) |>
#'   # fit the model:
#'   fit()
#' @export
addLsp <- function(mixedPenalty,
                   regularized,
                   lambdas,
                   thetas){
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  if(!is(mixedPenalty, "mixedPenalty"))
    stop("mixedPenalty must be of class mixedPenalty. ",
         "These models can be created with the regularize() function.")
  
  tps <- expand.grid(lambda = lambdas,
                     theta = thetas,
                     alpha = 0)
  
  mixedPenalty$penalties <- c(
    mixedPenalty$penalties,
    list(list(regularized = regularized, 
              tps = tps,
              penaltyType = .penaltyTypes("lsp"))
    )
  )
  
  return(mixedPenalty)
}

#' addMcp
#' 
#' Implements mcp regularization for structural equation models.
#' The penalty function is given by:
#' \ifelse{html}{\deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| - x_j^2/(2\theta) & \text{if } |x_j| \leq \theta\lambda\\
#' \theta\lambda^2/2 & \text{if } |x_j| > \lambda\theta
#' \end{cases}} where \eqn{\theta > 0}.}{
#' Equation Omitted in Pdf Documentation.}
#' 
#' Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
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
#' @param mixedPenalty model of class mixedPenalty created with the mixedPenalty function (see ?mixedPenalty) 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @returns Model of class mixedPenalty. Use the fit() - function to fit the model
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
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addMcp(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,.1),
#'           thetas = 2.3) |>
#'   # fit the model:
#'   fit()
#' @export
addMcp <- function(mixedPenalty,
                   regularized,
                   lambdas,
                   thetas){
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  if(!is(mixedPenalty, "mixedPenalty"))
    stop("mixedPenalty must be of class mixedPenalty. ",
         "These models can be created with the regularize() function.")
  
  tps <- expand.grid(lambda = lambdas,
                     theta = thetas,
                     alpha = 0)
  
  mixedPenalty$penalties <- c(
    mixedPenalty$penalties,
    list(list(regularized = regularized, 
              tps = tps,
              penaltyType = .penaltyTypes("mcp"))
    )
  )
  
  return(mixedPenalty)
}

#' addScad
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
#' @param mixedPenalty model of class mixedPenalty created with the mixedPenalty function (see ?mixedPenalty) 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @returns Model of class mixedPenalty. Use the fit() - function to fit the model
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
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addScad(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,.1),
#'           thetas = 3.1) |>
#'   # fit the model:
#'   fit()
#' @export
addScad <- function(mixedPenalty,
                    regularized,
                    lambdas,
                    thetas){
  
  if(any(thetas <= 2)) stop("Theta must be > 2")
  
  if(!is(mixedPenalty, "mixedPenalty"))
    stop("mixedPenalty must be of class mixedPenalty. ",
         "These models can be created with the regularize() function.")
  
  tps <- expand.grid(lambda = lambdas,
                     theta = thetas,
                     alpha = 0)
  
  mixedPenalty$penalties <- c(
    mixedPenalty$penalties,
    list(list(regularized = regularized, 
              tps = tps,
              penaltyType = .penaltyTypes("scad"))
    )
  )
  
  return(mixedPenalty)
}


#' addElasticNet
#' 
#' Adds an elastic net penalty to specified parameters.
#' The penalty function is given by:
#' \deqn{p( x_j) = \alpha\lambda|x_j| + (1-\alpha)\lambda x_j^2}
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
#' @param mixedPenalty model of class mixedPenalty created with the mixedPenalty function (see ?mixedPenalty) 
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param alphas numeric vector: values for the tuning parameter alpha. Set to 1 for lasso and
#' to zero for ridge. Anything in between is an elastic net penalty.
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param weights can be used to give different weights to the different parameters
#' @returns Model of class mixedPenalty. Use the fit() - function to fit the model
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
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addElasticNet(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,.1),
#'           alphas = .4) |>
#'   # fit the model:
#'   fit()
#' @export
addElasticNet <- function(mixedPenalty,
                          regularized,
                          alphas,
                          lambdas,
                          weights = 1){
  
  if(!is(mixedPenalty, "mixedPenalty"))
    stop("mixedPenalty must be of class mixedPenalty. ",
         "These models can be created with the regularize() function.")
  if(!mixedPenalty$method == "glmnet"){
    notes <- c(notes,
               "Mixed penalty with addElasticNet is only supported for method = 'glmnet'. Switching optimizer to glmnet")
    mixedPenalty$method = "glmnet"
    mixedPenalty$control = controlGlmnet()
  }
  
  if(any(alphas > 1) || any(alphas < 0))
    stop("alphas must be in the interval [0,1]")
  
  tps <- expand.grid(alpha = alphas,
                     lambda = lambdas)
  
  if((length(weights) != 1) && (length(weights) != length(regularized)))
    stop("Weights argument must be of length 1 or of length(regularized).")
  
  mixedPenalty$penalties <- c(
    mixedPenalty$penalties,
    list(list(regularized = regularized,
              weight = weights,
              tps = tps,
              penaltyType = "elasticNet"))
  )
  
  return(mixedPenalty)
}

#' fit
#' 
#' Optimizes an object with mixed penalty. See ?mixedPenalty for more details.
#' 
#' @param mixedPenalty object of class mixedPenalty. This object can be created
#' with the mixedPenalty function. Penalties can be added with the addCappedL1, addElastiNet,
#' addLasso, addLsp, addMcp, and addScad functions.
#' @return throws error in case of undefined penalty combinations.
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
#' # We can add mixed penalties as follows:
#' 
#' regularized <- lavaanModel |>
#'   # create template for regularized model with mixed penalty:
#'   mixedPenalty() |>
#'   # add penalty on loadings l6 - l10:
#'   addElasticNet(regularized = paste0("l", 11:15), 
#'           lambdas = seq(0,1,.1),
#'           alphas = .4) |>
#'   # fit the model:
#'   fit()
#' @export 
fit <- function(mixedPenalty){
  
  .checkPenalties(mixedPenalty)
  
  if(.useElasticNet(mixedPenalty))
    return(.fitElasticNetMix(mixedPenalty))
  
  return(.fitMix(mixedPenalty))
  
  stop("Unknown method used.")
}

#' .checkPenalties
#' 
#' Internal function to check a mixedPenalty object
#' @param mixedPenalty object of class mixedPenalty. This object can be created
#' with the mixedPenalty function. Penalties can be added with the addCappedL1,
#' addLasso, addLsp, addMcp, and addScad functions.
.checkPenalties <- function(mixedPenalty){
  
  if(!is(mixedPenalty, "mixedPenalty")){
    stop("mixedPenalty must be of class mixedPenalty.")
  }
  
  if(length(mixedPenalty$penalties) == 0)
    stop("Penaties are missing")
  
}

#' .useElasticNet
#' 
#' Internal function checking if elastic net is used
#' @param mixedPenalty object of class mixedPenalty. This object can be created
#' with the mixedPenalty function. Penalties can be added with the addCappedL1,
#' addLasso, addLsp, addMcp, and addScad functions.
#' @return TRUE if elastic net, FALSE otherwise
.useElasticNet <- function(mixedPenalty){
  hasElasticNet <- FALSE
  hasNonElasticNet <- FALSE
  for(i in 1:length(mixedPenalty$penalties)){
    if(mixedPenalty$penalties[[i]]$penaltyType == "elasticNet")
      hasElasticNet <- TRUE
    if(mixedPenalty$penalties[[i]]$penaltyType != "elasticNet")
      hasNonElasticNet <- TRUE
  }
  
  if(hasElasticNet && mixedPenalty$method == "ista")
    stop("Ista currently does not support elasticNet type penalties.")
  if(hasElasticNet && hasNonElasticNet)
    stop("lessSEM cannot mix penalties of type elastic net and other penalties.")
  
  return(hasElasticNet)
}

#' .fitMix
#' 
#' Optimizes an object with mixed penalty. See ?mixedPenalty for more details.
#' 
#' @param mixedPenalty object of class mixedPenalty. This object can be created
#' with the mixedPenalty function. Penalties can be added with the addCappedL1, addElastiNet,
#' addLasso, addLsp, addMcp, and addScad functions.
#' @return object of class regularizedSEMMixedPenalty
#' @keywords internal 
.fitMix <- function(mixedPenalty){
  
  notes <- c("Notes:")
  
  lavaanModel = mixedPenalty$lavaanModel
  method = mixedPenalty$method
  modifyModel = mixedPenalty$modifyModel
  control = mixedPenalty$control
  
  inputArguments <- as.list(environment())
  
  if(method == "ista" && !is(control, "controlIsta")) 
    stop("control must be of class controlIsta. See ?controlIsta.")
  
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet")
  
  .checkLavaanModel(lavaanModel = lavaanModel)
  
  ### initialize model ####
  startingValues <- control$startingValues
  
  ### initialize model ####
  if(is(lavaanModel, "lavaan")){
    SEM <- .initializeSEMForRegularization(lavaanModel = lavaanModel,
                                           startingValues = startingValues,
                                           modifyModel = modifyModel)
  }else if(is.vector(lavaanModel)){
    SEM <- .initializeMultiGroupSEMForRegularization(lavaanModels = lavaanModel,
                                                     startingValues = startingValues,
                                                     modifyModel = modifyModel)
    notes <- c(notes,
    "Multi-group model. Switching initialHessian from 'lavaan' to 'compute'"
    )
    control$initialHessian <- "compute"
  }
  
  # check if we have a likelihood objective function:
  usesLikelihood <- all(SEM$getEstimator() == "fiml")
  
  N <- SEM$sampleSize
  
  # get parameters in raw form
  startingValues <- .getParameters(SEM, raw = TRUE)
  rawParameters <- .getParameters(SEM, raw = TRUE)
  
  # set weights
  weights <- rep(0, length(startingValues))
  names(weights) <- names(startingValues)
  penaltyType <- rep(0, length(startingValues))
  names(penaltyType) <- names(startingValues)
  
  tpRows <- vector("list", length(mixedPenalty$penalties))
  names(tpRows) <- paste0("penalty", 1:length(tpRows))
  
  for(p in 1:length(mixedPenalty$penalties)){
    if(any(weights[mixedPenalty$penalties[[p]]$regularized] != 0))
      stop("It seems like some parameters appeared in multiple penalties. This is not supported. ",
           "The following parameters were found in multiple penalties: ",
           paste0(mixedPenalty$penalties[[p]]$regularized[weights[mixedPenalty$penalties[[p]]$regularized] != 0], sep = ","))
    
    if(all(mixedPenalty$penalties[[p]]$penalty == .penaltyTypes("lasso"))){
      weights[mixedPenalty$penalties[[p]]$regularized] <- mixedPenalty$penalties[[p]]$weight
    }else{
      weights[mixedPenalty$penalties[[p]]$regularized] <- 1
    }
    
    penaltyType[mixedPenalty$penalties[[p]]$regularized] <- mixedPenalty$penalties[[p]]$penaltyType
    
    tpRows[[p]] <- 1:nrow(mixedPenalty$penalties[[p]]$tp)
  }
  
  tpGrid <- unique(expand.grid(tpRows))
  
  lambda <- theta <- alpha <- matrix(0, 
                                     nrow = nrow(tpGrid), 
                                     ncol = length(startingValues),
                                     dimnames = list(NULL, names(startingValues)))
  for(p in 1:ncol(tpGrid)){
    
    lambda[,mixedPenalty$penalties[[p]]$regularized] <-  matrix(
      rep(mixedPenalty$penalties[[p]]$tps$lambda[tpGrid[,p]], 
          length(mixedPenalty$penalties[[p]]$regularized)),
      ncol = length(mixedPenalty$penalties[[p]]$regularized))
    
    
    theta[,mixedPenalty$penalties[[p]]$regularized] <-  matrix(
      rep(mixedPenalty$penalties[[p]]$tps$theta[tpGrid[,p]], 
          length(mixedPenalty$penalties[[p]]$regularized)),
      ncol = length(mixedPenalty$penalties[[p]]$regularized))
    
    alpha[,mixedPenalty$penalties[[p]]$regularized] <-  matrix(
      rep(mixedPenalty$penalties[[p]]$tps$alpha[tpGrid[,p]], 
          length(mixedPenalty$penalties[[p]]$regularized)),
      ncol = length(mixedPenalty$penalties[[p]]$regularized))
    
  }
  
  #### prepare regularized model object ####
  if(method == "ista"){
    controlIntern <- list(
      L0 = control$L0,
      eta = control$eta,
      accelerate = control$accelerate,
      maxIterOut = control$maxIterOut,
      maxIterIn = control$maxIterIn,
      breakOuter = control$breakOuter,
      convCritInner = control$convCritInner,
      sigma = control$sigma,
      stepSizeInheritance = control$stepSizeInheritance,
      verbose = control$verbose
    )
    
    
    
    if(is(SEM, "Rcpp_SEMCpp")){
      regularizedModel <- new(istaMixedPenaltySEM, 
                              weights, 
                              penaltyType,
                              controlIntern)
    }else if(is(SEM, "Rcpp_mgSEM")){
      regularizedModel <- new(istaMixedPenaltymgSEM, 
                              weights, 
                              penaltyType,
                              controlIntern)
    }
  }
  
  if(method == "glmnet"){
    #### glmnet requires an initial Hessian ####
    initialHessian <- .computeInitialHessian(initialHessian = control$initialHessian, 
                                            rawParameters = rawParameters, 
                                            lavaanModel = lavaanModel, 
                                            SEM = SEM,
                                            addMeans = modifyModel$addMeans,
                                            stepSize = control$stepSize,
                                            notes = notes)
    notes <- initialHessian$notes
    control$initialHessian <- initialHessian$initialHessian
    
    #### prepare regularized model object ####
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
    
    if(is(SEM, "Rcpp_SEMCpp")){
      regularizedModel <- new(glmnetMixedSEM, 
                              weights, 
                              penaltyType,
                              controlIntern)
    }else if(is(SEM, "Rcpp_mgSEM")){
      regularizedModel <- new(glmnetMixedMgSEM, 
                              weights, 
                              penaltyType,
                              controlIntern)
    }
  }
  
  #### define tuning parameters and prepare fit results ####
  
  fits <- data.frame(
    "tuningParameterConfiguration" = 1:nrow(tpGrid),
    "objectiveValue" = NA,
    "regObjectiveValue" = NA,
    "m2LL" = NA,
    "regM2LL"= NA,
    "nonZeroParameters" = NA,
    "convergence" = NA
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,
                                             nrow = nrow(tpGrid), 
                                             ncol = length(rawParameters)))
  colnames(parameterEstimates) <- names(rawParameters)
  parameterEstimates <- cbind(
    "tuningParameterConfiguration" = 1:nrow(tpGrid),
    parameterEstimates
  )
  
  if(!is.null(modifyModel$transformations)){
    transformationsAre <- .getParameters(SEM,
                                         raw = FALSE,
                                         transformations = TRUE)
    transformationsAre <- transformationsAre[!names(transformationsAre)%in%names(rawParameters)]
    
    transformations <- as.data.frame(matrix(NA,
                                            nrow = nrow(tpGrid), 
                                            ncol = length(transformationsAre)))
    colnames(transformations) <- names(transformationsAre)
    
    transformations <- cbind(
      "tuningParameterConfiguration" = 1:nrow(tpGrid),
      transformations
    )
  }else{
    transformations <- data.frame()
  }
  
  Hessians <- list(NULL)
  
  # save model implied matrices
  implied <- list(means = vector("list", nrow(tpGrid)),
                  covariances = vector("list", nrow(tpGrid)))
  
  #### print progress ####
  if(control$verbose == 0){
    progressbar = utils::txtProgressBar(min = 0, 
                                        max = nrow(tpGrid), 
                                        initial = 0, 
                                        style = 3)
  }
  
  #### Iterate over all tuning parameter combinations and fit models ####
  
  for(it in 1:nrow(tpGrid)){
    if(control$verbose == 0){
      utils::setTxtProgressBar(progressbar,it)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tpGrid),"]\n"))
    }
    
    result <- try(regularizedModel$optimize(rawParameters,
                                            SEM,
                                            lambda[it,],
                                            theta[it,],
                                            alpha[it,])
    )
    
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(rawParameters[weights[names(rawParameters)] != 0] == 0)
    
    fits$regObjectiveValue[it] <- result$fit
    if(usesLikelihood)
      fits$regM2LL[it] <- fits$regObjectiveValue[it]
    
    fits$convergence[it] <- result$convergence
    
    # get unregularized fit:
    SEM <- .setParameters(SEM, 
                          names(rawParameters), 
                          values = rawParameters, 
                          raw = TRUE)
    
    fits$objectiveValue[it] <- SEM$fit()
    if(usesLikelihood)
      fits$m2LL[it] <- fits$objectiveValue[it]
    
    # transform internal parameter representation to "natural" form
    transformedParameters <- .getParameters(SEM,
                                            raw = FALSE)
    parameterEstimates[it, 
                       names(rawParameters)] <- transformedParameters[names(rawParameters)]
    
    if(!is.null(modifyModel$transformations)){
      transformationsAre <- .getParameters(SEM,
                                           raw = FALSE,
                                           transformations = TRUE)
      transformationsAre <- transformationsAre[!names(transformationsAre)%in%names(rawParameters)]
      transformations[it,
                      names(transformationsAre)] <- transformationsAre
    }
    
    # save implied
    if(is(SEM, "Rcpp_SEMCpp") & control$saveDetails){
      implied$means[[it]] <- SEM$impliedMeans
      implied$covariances[[it]] <- SEM$impliedCovariance
      
      rownames(implied$means[[it]]) <- SEM$manifestNames
      dimnames(implied$covariances[[it]]) <- list(SEM$manifestNames,
                                                  SEM$manifestNames)
    }else if(is(SEM, "Rcpp_mgSEM")){
      implied$means[[it]] <- NULL
      implied$covariances[[it]] <- NULL
    }
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning(paste0("Fit for tuningParameterConfiguration = ", it,
                     " resulted in Error!"))
      
      if(is(lavaanModel, "lavaan")){
        SEM <- .initializeSEMForRegularization(lavaanModel = lavaanModel,
                                               startingValues = startingValues,
                                               modifyModel = modifyModel)
      }else if(is.vector(lavaanModel)){
        SEM <- .initializeMultiGroupSEMForRegularization(lavaanModels = lavaanModel,
                                                         startingValues = startingValues,
                                                         modifyModel = modifyModel)
      }
      
    }else{
      
      if(method == "glmnet"){
        # set Hessian for next iteration
        regularizedModel$setHessian(result$Hessian)
      }
      
    }
    
  }
  
  if(is(SEM, "Rcpp_SEMCpp")) internalOptimization <- list(
    "implied" = implied,
    "HessiansOfDifferentiablePart" = Hessians,
    "functionCalls" = SEM$functionCalls,
    "gradientCalls" = SEM$gradientCalls,
    "N" = SEM$sampleSize
  )
  if(is(SEM, "Rcpp_mgSEM")) internalOptimization <- list(
    "implied" = implied,
    "HessiansOfDifferentiablePart" = Hessians,
    "functionCalls" = NA,
    "gradientCalls" = NA,
    "N" = SEM$sampleSize
  )
  
  notes <- unique(notes)
  if(length(notes) > 1){
    cat("\n")
    rlang::inform(
      notes
    )
  }
  
  results <- new("regularizedSEMMixedPenalty",
                 penalty =  sapply(penaltyType, .penaltyTypes),
                 tuningParameterConfigurations = list(
                   lambda = cbind("configuration" = 1:nrow(tpGrid),
                                  lambda),
                   theta = cbind("configuration" = 1:nrow(tpGrid),
                                 theta),
                   alpha = cbind("configuration" = 1:nrow(tpGrid),
                                 alpha)
                 ),
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 transformations = transformations,
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments,
                 notes = notes)
  
  return(results)
  
}

#' .fitElasticNetMix
#' 
#' Optimizes an object with mixed penalty. See ?mixedPenalty for more details.
#' 
#' @param mixedPenalty object of class mixedPenalty. This object can be created
#' with the mixedPenalty function. Penalties can be added with the addCappedL1, addElastiNet,
#' addLasso, addLsp, addMcp, and addScad functions.
#' @return object of class regularizedSEMMixedPenalty
#' @keywords internal 
.fitElasticNetMix <- function(mixedPenalty){
  
  notes <- c("Notes:")
  
  lavaanModel = mixedPenalty$lavaanModel
  method = mixedPenalty$method
  modifyModel = mixedPenalty$modifyModel
  control = mixedPenalty$control
  
  inputArguments <- as.list(environment())
  
  if(! method %in% c("glmnet")) 
    stop("Currently ony methods = 'glmnet' is supported for mixtures of elastic net penalties.")
  
  if(method == "glmnet" && !is(control, "controlGlmnet")) 
    stop("control must be of class controlGlmnet See ?controlGlmnet.")
  
  .checkLavaanModel(lavaanModel = lavaanModel)
  
  ### initialize model ####
  startingValues <- control$startingValues
  
  ### initialize model ####
  if(is(lavaanModel, "lavaan")){
    SEM <- .initializeSEMForRegularization(lavaanModel = lavaanModel,
                                           startingValues = startingValues,
                                           modifyModel = modifyModel)
  }else if(is.vector(lavaanModel)){
    SEM <- .initializeMultiGroupSEMForRegularization(lavaanModels = lavaanModel,
                                                     startingValues = startingValues,
                                                     modifyModel = modifyModel)
  }
  
  N <- SEM$sampleSize
  
  # get parameters in raw form
  startingValues <- .getParameters(SEM, raw = TRUE)
  rawParameters <- .getParameters(SEM, raw = TRUE)
  
  # set weights
  weights <- rep(0, length(startingValues))
  names(weights) <- names(startingValues)
  penaltyType <- rep("none", length(startingValues))
  names(penaltyType) <- names(startingValues)
  
  tpRows <- vector("list", length(mixedPenalty$penalties))
  names(tpRows) <- paste0("penalty", 1:length(tpRows))
  
  for(p in 1:length(mixedPenalty$penalties)){
    if(any(weights[mixedPenalty$penalties[[p]]$regularized] != 0))
      stop("It seems like some parameters appeared in multiple penalties. This is not supported. ",
           "The following parameters were found in multiple penalties: ",
           paste0(mixedPenalty$penalties[[p]]$regularized[weights[mixedPenalty$penalties[[p]]$regularized] != 0], sep = ","))
    
    weights[mixedPenalty$penalties[[p]]$regularized] <- mixedPenalty$penalties[[p]]$weight
    
    penaltyType[mixedPenalty$penalties[[p]]$regularized] <- "elasticNet"
    
    tpRows[[p]] <- 1:nrow(mixedPenalty$penalties[[p]]$tp)
  }
  
  tpGrid <- unique(expand.grid(tpRows))
  
  lambda <- alpha <- matrix(0, 
                            nrow = nrow(tpGrid), 
                            ncol = length(startingValues),
                            dimnames = list(NULL, names(startingValues)))
  for(p in 1:ncol(tpGrid)){
    
    lambda[,mixedPenalty$penalties[[p]]$regularized] <-  matrix(
      rep(mixedPenalty$penalties[[p]]$tps$lambda[tpGrid[,p]], 
          length(mixedPenalty$penalties[[p]]$regularized)),
      ncol = length(mixedPenalty$penalties[[p]]$regularized))
    
    alpha[,mixedPenalty$penalties[[p]]$regularized] <-  matrix(
      rep(mixedPenalty$penalties[[p]]$tps$alpha[tpGrid[,p]], 
          length(mixedPenalty$penalties[[p]]$regularized)),
      ncol = length(mixedPenalty$penalties[[p]]$regularized))
    
  }
  
  #### glmnet requires an initial Hessian ####
  initialHessian <- .computeInitialHessian(initialHessian = control$initialHessian, 
                                          rawParameters = rawParameters, 
                                          lavaanModel = lavaanModel, 
                                          SEM = SEM,
                                          addMeans = modifyModel$addMeans,
                                          stepSize = control$stepSize,
                                          notes = notes)
  notes <- initialHessian$notes
  control$initialHessian <- initialHessian$initialHessian
  
  #### prepare regularized model object ####
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
  
  if(is(SEM, "Rcpp_SEMCpp")){
    regularizedModel <- new(glmnetEnetSEM, 
                            weights, 
                            controlIntern)
  }else if(is(SEM, "Rcpp_mgSEM")){
    regularizedModel <- new(glmnetEnetMgSEM, 
                            weights, 
                            controlIntern)
  }
  
  
  #### define tuning parameters and prepare fit results ####
  
  fits <- data.frame(
    "tuningParameterConfiguration" = 1:nrow(tpGrid),
    "m2LL" = NA,
    "regM2LL"= NA,
    "nonZeroParameters" = NA,
    "convergence" = NA
  )
  
  parameterEstimates <- as.data.frame(matrix(NA,
                                             nrow = nrow(tpGrid), 
                                             ncol = length(rawParameters)))
  colnames(parameterEstimates) <- names(rawParameters)
  parameterEstimates <- cbind(
    "tuningParameterConfiguration" = 1:nrow(tpGrid),
    parameterEstimates
  )
  
  if(!is.null(modifyModel$transformations)){
    transformationsAre <- .getParameters(SEM,
                                         raw = FALSE,
                                         transformations = TRUE)
    transformationsAre <- transformationsAre[!names(transformationsAre)%in%names(rawParameters)]
    
    transformations <- as.data.frame(matrix(NA,
                                            nrow = nrow(tpGrid), 
                                            ncol = length(transformationsAre)))
    colnames(transformations) <- names(transformationsAre)
    
    transformations <- cbind(
      "tuningParameterConfiguration" = 1:nrow(tpGrid),
      transformations
    )
  }else{
    transformations <- data.frame()
  }
  
  if(control$saveDetails){
    Hessians <- list(
      "lambda" = lambda,
      "alpha" = alpha,
      "Hessian" = lapply(1:nrow(tpGrid), 
                         matrix, 
                         data= NA, 
                         nrow=nrow(control$initialHessian), 
                         ncol=ncol(control$initialHessian))
    )
  }else{
    Hessians <- list(NULL)
  }
  
  # save model implied matrices
  implied <- list(means = vector("list", nrow(tpGrid)),
                  covariances = vector("list", nrow(tpGrid)))
  
  #### print progress ####
  if(control$verbose == 0){
    progressbar = utils::txtProgressBar(min = 0, 
                                        max = nrow(tpGrid), 
                                        initial = 0, 
                                        style = 3)
  }
  
  #### Iterate over all tuning parameter combinations and fit models ####
  
  for(it in 1:nrow(tpGrid)){
    if(control$verbose == 0){
      utils::setTxtProgressBar(progressbar,it)
    }else{
      cat(paste0("\nIteration [", it, "/", nrow(tpGrid),"]\n"))
    }
    
    result <- try(regularizedModel$optimize(rawParameters,
                                            SEM,
                                            lambda[it,],
                                            alpha[it,]
    )
    )
    
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(rawParameters[weights[names(rawParameters)] != 0] == 0)
    fits$regM2LL[it] <- result$fit
    fits$convergence[it] <- result$convergence
    
    # get unregularized fit:
    SEM <- .setParameters(SEM, 
                          names(rawParameters), 
                          values = rawParameters, 
                          raw = TRUE)
    fits$m2LL[it] <- SEM$fit()
    # transform internal parameter representation to "natural" form
    transformedParameters <- .getParameters(SEM,
                                            raw = FALSE)
    parameterEstimates[it, 
                       names(rawParameters)] <- transformedParameters[names(rawParameters)]
    
    if(!is.null(modifyModel$transformations)){
      transformationsAre <- .getParameters(SEM,
                                           raw = FALSE,
                                           transformations = TRUE)
      transformationsAre <- transformationsAre[!names(transformationsAre)%in%names(rawParameters)]
      transformations[it,
                      names(transformationsAre)] <- transformationsAre
    }
    
    # save implied
    if(is(SEM, "Rcpp_SEMCpp") & control$saveDetails){
      implied$means[[it]] <- SEM$impliedMeans
      implied$covariances[[it]] <- SEM$impliedCovariance
      
      rownames(implied$means[[it]]) <- SEM$manifestNames
      dimnames(implied$covariances[[it]]) <- list(SEM$manifestNames,
                                                  SEM$manifestNames)
    }else if(is(SEM, "Rcpp_mgSEM")){
      implied$means[[it]] <- NULL
      implied$covariances[[it]] <- NULL
    }
    
    # set initial values for next iteration
    if(is(SEM, "try-Error")){
      # reset
      warning(paste0("Fit for tuningParameterConfiguration = ", it,
                     " resulted in Error!"))
      
      if(is(lavaanModel, "lavaan")){
        SEM <- .initializeSEMForRegularization(lavaanModel = lavaanModel,
                                               startingValues = startingValues,
                                               modifyModel = modifyModel)
      }else if(is.vector(lavaanModel)){
        SEM <- .initializeMultiGroupSEMForRegularization(lavaanModels = lavaanModel,
                                                         startingValues = startingValues,
                                                         modifyModel = modifyModel)
      }
      
    }else{
      
      if(control$saveDetails) Hessians$Hessian[[it]] <- result$Hessian
      
      # set Hessian for next iteration
      regularizedModel$setHessian(result$Hessian)
      
    }
    
  }
  
  if(is(SEM, "Rcpp_SEMCpp")) internalOptimization <- list(
    "implied" = implied,
    "HessiansOfDifferentiablePart" = Hessians,
    "functionCalls" = SEM$functionCalls,
    "gradientCalls" = SEM$gradientCalls,
    "N" = SEM$sampleSize
  )
  if(is(SEM, "Rcpp_mgSEM")) internalOptimization <- list(
    "implied" = implied,
    "HessiansOfDifferentiablePart" = Hessians,
    "functionCalls" = NA,
    "gradientCalls" = NA,
    "N" = SEM$sampleSize
  )
  
  notes <- unique(notes)
  if(length(notes) > 1){
    cat("\n")
    rlang::inform(
      notes
    )
  }
  
  results <- new("regularizedSEMMixedPenalty",
                 penalty =  penaltyType,
                 tuningParameterConfigurations = list(
                   lambda = cbind("configuration" = 1:nrow(tpGrid),
                                  lambda),
                   alpha = cbind("configuration" = 1:nrow(tpGrid),
                                 alpha)
                 ),
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 transformations = transformations,
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments,
                 notes = notes)
  
  return(results)
  
}

#' .penaltyTypes
#' 
#' translates the penalty from a numeric value to the character or from
#' the character to the numeric value. The numeric value is used by the C++ backend.
#' @param penalty either a number or the name of the penalty
#' @return number corresponding to one of the penalties
#' @keywords internal
.penaltyTypes <- function(penalty){
  
  if(is.character(penalty)){
    
    if(penalty == "none"){
      return(0)
    }
    
    if(penalty == "cappedL1"){
      return(1)
    }
    
    if(penalty == "lasso"){
      return(2)
    }
    
    if(penalty == "lsp"){
      return(3)
    }
    
    if(penalty == "mcp"){
      return(4)
    }
    
    if(penalty == "scad"){
      return(5)
    }
    
    stop("Unknown penalty.")
  }
  
  if(is.numeric(penalty)){
    
    if(penalty == 0){
      return("none")
    }
    
    if(penalty == 1){
      return("cappedL1")
    }
    
    if(penalty == 2){
      return("lasso")
    }
    
    if(penalty == 3){
      return("lsp")
    }
    
    if(penalty == 4){
      return("mcp")
    }
    
    if(penalty == 5){
      return("scad")
    }
    
    stop("Unknown penalty.")
  }
  
  stop("Unknown penalty.")
}