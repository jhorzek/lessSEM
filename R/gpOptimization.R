#' gpLasso
#' 
#' Optimizers in lessSEM can also be used by other packages. 
#' gpLasso allows for combining lasso regularization and custom fitting functions
#' with an interface similar to optim.  
#' 
#' For more details on GLMNET, see:
#' 
#' 1. Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' 2. Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
#' 3. Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' 1. Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' 2. Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#'  
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr R function which takes the parameters as input and returns the 
#' gradients of the objective function. If set to NULL, numDeriv will be used
#' to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param additionalArguments additional argument passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @md
#' @examples 
#' # This example shows how to use the optimizers
#' # for other objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(lessSEM)
#' set.seed(123)
#' 
#' # first, we simulate data for our
#' # linear regression. The model is given by
#' # y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 
#' n <- 100
#' df <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   x3 = rnorm(n),
#'   x4 = rnorm(n),
#'   x5 = rnorm(n)
#' )
#' df$y <- 0*df$x1 + .2*df$x2 + .3*df$x3 + .4*df$x4 + .5*df$x5 + rnorm(n,0,.3) 
#' 
#' # First, we must construct a fiting function
#' # which returns a single value. We will use
#' # the residual sum squared as fitting function.
#' # The optimizer expects that this function takes
#' # EXACTLY three arguments:
#' # 1) a vector with parameter values (par)
#' # 2) a vector with labels for these parameters (parameterLabels)
#' # 3) an additional argument which allows for passing anything
#' # else your function needs to run to the function.
#' 
#' # Let's start setting up the fitting function:
#' fittingFunction <- function(par, parameterLabels, additionalArguments){
#'   # Our function needs the observed data 
#'   # y and X. These are not part of the first two arguments
#'   # and must therefore be passed in the function with the 
#'   # third argument (additionalArguments). Here, we will use a list
#'   # and pass the arguments in through that list:
#'   pred <- additionalArguments$X %*% matrix(par, ncol = 1)
#'   sse <- sum((additionalArguments$y - pred)^2)
#'   return(sse)
#' }
#' # Now we need to construct the list additionalArguments used
#' # by fittingFunction
#' additionalArguments <- list(X = as.matrix(df[,grepl("x", colnames(df))]), 
#'                             y = df$y)
#' 
#' 
#' # let's define the starting values:
#' b <- rep(0,5)
#' names(b) <- paste0("b",1:5)
#' # names of regularized parameters
#' regularized <- paste0("b",1:5)
#' 
#' # optimize
#' lassoPen <- gpLasso(
#'   par = b, 
#'   regularized = regularized, 
#'   fn = fittingFunction, 
#'   nLambdas = 100, 
#'   additionalArguments = additionalArguments
#' )
#' plot(lassoPen)
#' AIC(lassoPen)
#' @export
gpLasso <- function(par,
                    regularized,
                    fn,
                    gr = NULL,
                    lambdas = NULL,
                    nLambdas = NULL,
                    additionalArguments,
                    method = "ista", 
                    control = controlIsta()
){
  
  weights <- par
  weights[] <- 0
  weights[regularized] <- 1
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  result <- gpElasticNet(
    par = par,
    fn = fn,
    gr = gr,
    weights = weights,
    lambdas = lambdas,
    nLambdas = nLambdas,
    alphas = 1,
    additionalArguments = additionalArguments,
    method = method,
    control = control
  )
  return(result)
  
}

#' gpAdaptiveLasso
#' 
#' Optimizers in lessSEM can also be used by other packages. 
#' gpAdaptiveLasso allows for combining adaptive lasso regularization and custom fitting functions
#' with an interface similar to optim.  
#' 
#' For more details on GLMNET, see:
#' 
#' 1. Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' 2. Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
#' 3. Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' 1. Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' 2. Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#'  
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr R function which takes the parameters as input and returns the 
#' gradients of the objective function. If set to NULL, numDeriv will be used
#' to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param weights labeled vector with adaptive lasso weights. NULL will use 1/abs(par)
#' @param additionalArguments additional argument passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @md
#' @examples 
#' # This example shows how to use the optimizers
#' # for other objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(lessSEM)
#' set.seed(123)
#' 
#' # first, we simulate data for our
#' # linear regression. The model is given by
#' # y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 
#' n <- 100
#' df <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   x3 = rnorm(n),
#'   x4 = rnorm(n),
#'   x5 = rnorm(n)
#' )
#' df$y <- 0*df$x1 + .2*df$x2 + .3*df$x3 + .4*df$x4 + .5*df$x5 + rnorm(n,0,.3) 
#' 
#' # First, we must construct a fiting function
#' # which returns a single value. We will use
#' # the residual sum squared as fitting function.
#' # The optimizer expects that this function takes
#' # EXACTLY three arguments:
#' # 1) a vector with parameter values (par)
#' # 2) a vector with labels for these parameters (parameterLabels)
#' # 3) an additional argument which allows for passing anything
#' # else your function needs to run to the function.
#' 
#' # Let's start setting up the fitting function:
#' fittingFunction <- function(par, parameterLabels, additionalArguments){
#'   # Our function needs the observed data 
#'   # y and X. These are not part of the first two arguments
#'   # and must therefore be passed in the function with the 
#'   # third argument (additionalArguments). Here, we will use a list
#'   # and pass the arguments in through that list:
#'   pred <- additionalArguments$X %*% matrix(par, ncol = 1)
#'   sse <- sum((additionalArguments$y - pred)^2)
#'   return(sse)
#' }
#' # Now we need to construct the list additionalArguments used
#' # by fittingFunction
#' additionalArguments <- list(X = as.matrix(df[,grepl("x", colnames(df))]), 
#'                             y = df$y)
#' 
#' 
#' # let's define the starting values:
#' b <- rep(0,5)
#' names(b) <- paste0("b",1:5)
#' # names of regularized parameters
#' regularized <- paste0("b",1:5)
#' 
#' # for our weights, we need the unregularized parameter
#' # estimates
#' mle <- coef(lm(y~ 0+., data = df))
#' weights <- 1/abs(mle)
#' names(weights) <- paste0("b",1:5)
#' weights[!names(weights) %in% regularized] <- 0
#' 
#' # optimize
#' alasso <- gpAdaptiveLasso(
#'   par = b, 
#'   regularized = regularized, 
#'   weights = weights,
#'   fn = fittingFunction, 
#'   nLambdas = 100, 
#'   additionalArguments = additionalArguments
#' )
#' plot(alasso)
#' AIC(alasso)
#' @export
gpAdaptiveLasso <- function(par,
                            regularized,
                            weights = NULL,
                            fn,
                            gr = NULL,
                            lambdas = NULL,
                            nLambdas = NULL,
                            additionalArguments = NULL,
                            method = "ista", 
                            control = controlIsta()){
  if(is.null(weights)){
    weights <- 1/abs(par)
    weights[!names(weights) %in% regularized] <- 0
  }
  
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  result <- gpElasticNet(
    par = par,
    fn = fn,
    gr = gr,
    weights = weights,
    lambdas = lambdas,
    nLambdas = nLambdas,
    alphas = 1,
    additionalArguments = additionalArguments,
    method = method,
    control = control
  )
  return(result)
  
}

#' gpRidge
#' 
#' Optimizers in lessSEM can also be used by other packages. 
#' gpRidge allows for combining ridge regularization and custom fitting functions
#' with an interface similar to optim.  
#' 
#' For more details on GLMNET, see:
#' 
#' 1. Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' 2. Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
#' 3. Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' 1. Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' 2. Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#'  
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr R function which takes the parameters as input and returns the 
#' gradients of the objective function. If set to NULL, numDeriv will be used
#' to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param additionalArguments additional argument passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @md
#' @examples 
#' # This example shows how to use the optimizers
#' # for other objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(lessSEM)
#' set.seed(123)
#' 
#' # first, we simulate data for our
#' # linear regression. The model is given by
#' # y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 
#' n <- 100
#' df <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   x3 = rnorm(n),
#'   x4 = rnorm(n),
#'   x5 = rnorm(n)
#' )
#' df$y <- 0*df$x1 + .2*df$x2 + .3*df$x3 + .4*df$x4 + .5*df$x5 + rnorm(n,0,.3) 
#' 
#' # First, we must construct a fiting function
#' # which returns a single value. We will use
#' # the residual sum squared as fitting function.
#' # The optimizer expects that this function takes
#' # EXACTLY three arguments:
#' # 1) a vector with parameter values (par)
#' # 2) a vector with labels for these parameters (parameterLabels)
#' # 3) an additional argument which allows for passing anything
#' # else your function needs to run to the function.
#' 
#' # Let's start setting up the fitting function:
#' fittingFunction <- function(par, parameterLabels, additionalArguments){
#'   # Our function needs the observed data 
#'   # y and X. These are not part of the first two arguments
#'   # and must therefore be passed in the function with the 
#'   # third argument (additionalArguments). Here, we will use a list
#'   # and pass the arguments in through that list:
#'   pred <- additionalArguments$X %*% matrix(par, ncol = 1)
#'   sse <- sum((additionalArguments$y - pred)^2)
#'   return(sse)
#' }
#' # Now we need to construct the list additionalArguments used
#' # by fittingFunction
#' additionalArguments <- list(X = as.matrix(df[,grepl("x", colnames(df))]), 
#'                             y = df$y)
#' 
#' 
#' # let's define the starting values:
#' b <- rep(0,5)
#' names(b) <- paste0("b",1:5)
#' # names of regularized parameters
#' regularized <- paste0("b",1:5)
#' 
#' fittingFunction(b, names(b), additionalArguments)
#' 
#' # optimize
#' ridgeFit <- gpRidge(
#'   par = b, 
#'   regularized = regularized, 
#'   fn = fittingFunction, 
#'   lambdas = seq(0,100,1),
#'   additionalArguments = additionalArguments
#' )
#' plot(ridgeFit)
#' @export
gpRidge <- function(par,
                    regularized,
                    fn,
                    gr = NULL,
                    lambdas,
                    additionalArguments,
                    method = "ista", 
                    control = controlIsta()){
  
  weights <- par
  weights[] <- 0
  weights[regularized] <- 1
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  result <- gpElasticNet(
    par = par,
    fn = fn,
    gr = gr,
    weights = weights,
    lambdas = lambdas,
    nLambdas = NULL,
    alphas = 0,
    additionalArguments = additionalArguments,
    method = method,
    control = control
  )
  return(result)
  
}

#' gpElasticNet
#' 
#' Optimizers in lessSEM can also be used by other packages. 
#' gpElasticNet allows for combining elastic net regularization and custom 
#' fitting functions with an interface similar to optim.  
#' 
#' For more details on GLMNET, see:
#' 
#' 1. Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' 2. Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
#' 3. Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. 
#' The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421
#' 
#' For more details on ISTA, see:
#' 
#' 1. Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
#' 183–202. https://doi.org/10.1137/080716542
#' 2. Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
#' A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
#' Regularized Optimization Problems. Proceedings of the 30th International 
#' Conference on Machine Learning, 28(2)(2), 37–45.
#'  
#' @param par labeled vector with starting values
#' @param weights labeled vector with weights for each of the parameters in the 
#' model.
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param nLambdas alternative to lambda: If alpha = 1, lessSEM can automatically
#' compute the first lambda value which sets all regularized parameters to zero.
#' It will then generate nLambda values between 0 and the computed lambda.
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' in [0,1]. 0 = ridge, 1 = lasso.
#' @param additionalArguments additional argument passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related procedures
#' (currently gist).
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta() and controlGlmnet() functions.
#' @md
#' @examples
#' # This example shows how to use the optimizers
#' # for other objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
#' library(lessSEM)
#' set.seed(123)
#' 
#' # first, we simulate data for our
#' # linear regression. The model is given by
#' # y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 
#' n <- 100
#' df <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n),
#'   x3 = rnorm(n),
#'   x4 = rnorm(n),
#'   x5 = rnorm(n)
#' )
#' df$y <- 0*df$x1 + .2*df$x2 + .3*df$x3 + .4*df$x4 + .5*df$x5 + rnorm(n,0,.3) 
#' 
#' # First, we must construct a fiting function
#' # which returns a single value. We will use
#' # the residual sum squared as fitting function.
#' # The optimizer expects that this function takes
#' # EXACTLY three arguments:
#' # 1) a vector with parameter values (par)
#' # 2) a vector with labels for these parameters (parameterLabels)
#' # 3) an additional argument which allows for passing anything
#' # else your function needs to run to the function.
#' 
#' # Let's start setting up the fitting function:
#' fittingFunction <- function(par, parameterLabels, additionalArguments){
#'   # Our function needs the observed data 
#'   # y and X. These are not part of the first two arguments
#'   # and must therefore be passed in the function with the 
#'   # third argument (additionalArguments). Here, we will use a list
#'   # and pass the arguments in through that list:
#'   pred <- additionalArguments$X %*% matrix(par, ncol = 1)
#'   sse <- sum((additionalArguments$y - pred)^2)
#'   return(sse)
#' }
#' # Now we need to construct the list additionalArguments used
#' # by fittingFunction
#' additionalArguments <- list(X = as.matrix(df[,grepl("x", colnames(df))]), 
#'                             y = df$y)
#' 
#' 
#' # let's define the starting values:
#' b <- rep(0,5)
#' names(b) <- paste0("b",1:5)
#' # and the weights used by the elastic net
#' weights <- b
#' weights[] <- 1
#' 
#' # optimize
#' enet <- gpElasticNet(
#'   par = b, 
#'   weights = weights, 
#'   fn = fittingFunction, 
#'   lambdas = seq(0,100,1), 
#'   alphas = 1, 
#'   additionalArguments = additionalArguments
#' )
#' plot(enet)
#' AIC(enet)
#' @export
gpElasticNet <- function(par,
                         weights,
                         fn,
                         gr = NULL,
                         lambdas = NULL,
                         nLambdas = NULL,
                         alphas,
                         additionalArguments,
                         method = "ista", 
                         control = controlIsta()){
  
  inputArguments <- as.list(environment())
  
  parameterLabels <- names(par)
  if(is.null(parameterLabels)) stop("par has no labels.")
  rawParameters <- par
  
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
  
  # gradient function
  if(is.null(gr)){
    
    additionalArguments$fn <- fn
    
    gr <- function(par, parameterLabels, additionalArguments){
      grad <- numDeriv::grad(func = additionalArguments$fn, 
                             x = par, 
                             parameterLabels = parameterLabels,
                             additionalArguments = additionalArguments)
      names(grad) <- names(par)
      return(grad)
    }
    
  }
  
  # check functions
  check_fn <- fn(par, parameterLabels, additionalArguments)
  if(length(check_fn) != 1) stop("fn returns more than one element!")
  check_fn <- gr(par, parameterLabels, additionalArguments)
  if(length(check_fn) != length(par)) stop("gr has different length than par!")
  
  # make sure that the weights are in the correct order
  if(is.null(names(weights))) stop("weights must have the same names as the parameters")
  if(length(weights) != length(par)) stop("weights must be of the same length as the parameter vector.")
  if(any(!is.numeric(weights))) stop("weights must be numeric")
  weights <- weights[names(par)]
  
  #### glmnet requires an initial Hessian ####
  if(method == "glmnet"){
    initialHessian <- control$initialHessian
    if(is.matrix(initialHessian) && 
       nrow(initialHessian) == length(rawParameters) && 
       ncol(initialHessian) == length(rawParameters)){
      
    }else if(any(initialHessian == "compute")){
      
      initialHessian <- numDeriv::hessian(func = fn, x = par)
      
    }else if(length(initialHessian) == 1 && is.numeric(initialHessian)){
      initialHessian <- diag(initialHessian,length(rawParameters))
      rownames(initialHessian) <- names(rawParameters)
      colnames(initialHessian) <- names(rawParameters)
    }else{
      stop("Invalid initialHessian passed to GLMNET. See ?controlGLMNET for more information.")
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
      breakOuter = control$breakOuter,
      breakInner = control$breakInner,
      convergenceCriterion = control$convergenceCriterion, 
      verbose = control$verbose
    )
    
    regularizedModel <- new(glmnetEnetGeneralPurpose, 
                            weights, 
                            controlIntern)
    
  }else if(method == "ista"){
    
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
    
    regularizedModel <- new(istaEnetGeneralPurpose, 
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
    
    maxLambda <- gpGetMaxLambda(regularizedModel,
                                par,
                                fn,
                                gr,
                                additionalArguments,
                                weights)
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
    
    result <- try(
      regularizedModel$optimize( par,
                                 fn,
                                 gr,
                                 additionalArguments,
                                 lambda,
                                 alpha)
    )
    if(is(result, "try-error")) next
    
    rawParameters <- result$rawParameters
    parameterEstimates[it, names(rawParameters)] <- rawParameters[names(rawParameters)]
    
    fits$nonZeroParameters[it] <- length(rawParameters) - 
      sum(rawParameters[weights[names(rawParameters)] != 0] == 0)
    fits$regM2LL[it] <- result$fit
    fits$convergence[it] <- result$convergence
    fits$m2LL[it] <- fn(par, parameterLabels, additionalArguments)
    
    if(method == "glmnet" && control$saveHessian) 
      Hessians$Hessian[[it]] <- result$Hessian
    
    # set initial values for next iteration
    
    
    if(method == "glmnet"){
      if(control$saveHessian) Hessians$Hessian[[it]] <- result$Hessian
      
      # set Hessian for next iteration
      regularizedModel$setHessian(result$Hessian)
      
    }
    
  }
  
  internalOptimization <- list(
    "HessiansOfDifferentiablePart" = Hessians
  )
  
  results <- new("gpRegularized",
                 parameters = parameterEstimates,
                 fits = fits,
                 parameterLabels = names(rawParameters),
                 weights = weights,
                 regularized = names(weights)[weights!=0],
                 internalOptimization = internalOptimization,
                 inputArguments = inputArguments)
  
  return(results)
  
}