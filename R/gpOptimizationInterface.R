#' gpLasso 
#' 
#' Implements lasso regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda |x_j|}
#' Lasso regularization will set parameters to zero if \eqn{\lambda} is large enough
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with ...
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' Lasso regularization:
#' 
#' * Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical 
#' Society. Series B (Methodological), 58(1), 267–288.
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
#' # linear regression.
#' N <- 100 # number of persons
#' p <- 10 # number of predictors
#' X <- matrix(rnorm(N*p),	nrow = N, ncol = p) # design matrix
#' b <- c(rep(1,4), 
#'        rep(0,6)) # true regression weights
#' y <- X%*%matrix(b,ncol = 1) + rnorm(N,0,.2)
#' 
#' # First, we must construct a fiting function
#' # which returns a single value. We will use
#' # the residual sum squared as fitting function.
#' 
#' # Let's start setting up the fitting function:
#' fittingFunction <- function(par, y, X, N){
#'   # par is the parameter vector
#'   # y is the observed dependent variable
#'   # X is the design matrix
#'   # N is the sample size
#'   pred <- X %*% matrix(par, ncol = 1) #be explicit here: 
#'   # we need par to be a column vector
#'   sse <- sum((y - pred)^2)
#'   # we scale with .5/N to get the same results as glmnet
#'   return((.5/N)*sse)
#' }
#' 
#' # let's define the starting values:
#' b <- rep(0,p)
#' names(b) <- paste0("b", 1:length(b))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # optimize
#' lassoPen <- gpLasso(
#'   par = b, 
#'   regularized = regularized, 
#'   fn = fittingFunction, 
#'   nLambdas = 100, 
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' plot(lassoPen)
#' AIC(lassoPen)
#' 
#' # for comparison:
#' #library(glmnet)
#' #coef(glmnet(x = X,
#' #            y = y, 
#' #            lambda = lassoPen@fits$lambda[20],
#' #            intercept = FALSE,
#' #            standardize = FALSE))[,1]
#' #lassoPen@parameters[20,]
#' @export
gpLasso <- function(par,
                    regularized,
                    fn,
                    gr = NULL,
                    lambdas = NULL,
                    nLambdas = NULL,
                    ...,
                    method = "ista", 
                    control = controlIsta()
){
  
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1)
  }
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "lasso", 
                                            weights = regularized,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpAdaptiveLasso 
#' 
#' Implements adaptive lasso regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = p( x_j) = \frac{1}{w_j}\lambda| x_j|}
#' Adaptive lasso regularization will set parameters to zero if \eqn{\lambda} is large enough.
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take exactly 3 arguments:
#' 
#' 1. A vector with current parameter estimates
#' 2. A vector with labels of parameters
#' 3. Any additional arguments used by the function contained in a single list
#' (called additionalArguments below)
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' Adaptive lasso regularization:
#' 
#' * Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American Statistical Association, 
#' 101(476), 1418–1429. https://doi.org/10.1198/016214506000000735
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
                            ...,
                            method = "ista", 
                            control = controlIsta()){
  
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(is.null(weights)){
    weights <- 1/abs(par)
    weights[!names(weights) %in% regularized] <- 0
    message("Building weights based on par as weights = 1/abs(par).")
  }
  
  if(! all(regularized %in% names(weights))) stop(paste0(
    "You specified that the following parameters should be regularized:\n",
    paste0(regularized, collapse = ", "), 
    ". Not all of these parameters could be found in the model.\n",
    "The model has the following parameters:\n",
    names(weights)
  ))
  
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1)
  }
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "adaptiveLasso", 
                                            weights = regularized,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpRidge
#' 
#' Implements ridge regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda x_j^2}
#' Note that ridge regularization will not set any of the parameters to zero
#' but result in a shrinkage towards zero. 
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take exactly 3 arguments:
#' 
#' 1. A vector with current parameter estimates
#' 2. A vector with labels of parameters
#' 3. Any additional arguments used by the function contained in a single list
#' (called additionalArguments below)
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' Ridge regularization:
#' 
#' * Hoerl, A. E., & Kennard, R. W. (1970). Ridge Regression: Biased Estimation 
#' for Nonorthogonal Problems. Technometrics, 12(1), 55–67. 
#' https://doi.org/10.1080/00401706.1970.10488634
#' 
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
                    ...,
                    method = "ista", 
                    control = controlIsta()){
  
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 0)
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "ridge", 
                                            weights = regularized,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
  
}

#' gpElasticNet 
#' 
#' Implements elastic net regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = p( x_j) = \frac{1}{w_j}\lambda| x_j|}
#' Note that the elastic net combines ridge and lasso regularization. If \eqn{\alpha = 0}, 
#' the elastic net reduces to ridge regularization. If \eqn{\alpha = 1} it reduces
#' to lasso regularization. In between, elastic net is a compromise between the shrinkage of
#' the lasso and the ridge penalty. 
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take exactly 3 arguments:
#' 
#' 1. A vector with current parameter estimates
#' 2. A vector with labels of parameters
#' 3. Any additional arguments used by the function contained in a single list
#' (called additionalArguments below)
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' Elastic net regularization:
#' 
#' * Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. 
#' Journal of the Royal Statistical Society: Series B, 67(2), 301–320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
                         regularized,
                         fn,
                         gr = NULL,
                         lambdas,
                         alphas,
                         ...,
                         method = "ista", 
                         control = controlIsta()){
  
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = alphas)
  
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "elasticNet", 
                                            weights = regularized,
                                            tuningParameters = tuningParameters, 
                                            method = method,
                                            control = control
  )
  
  return(result)
}


#' gpCappedL1 
#' 
#' Implements cappedL1 regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \min(| x_j|, \theta)}
#' where \eqn{\theta > 0}. The cappedL1 penalty is identical to the lasso for 
#' parameters which are below \eqn{\theta} and identical to a constant for parameters
#' above \eqn{\theta}. As adding a constant to the fitting function will not change its
#' minimum, larger parameters can stay unregularized while smaller ones are set to zero.
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take exactly 3 arguments:
#' 
#' 1. A vector with current parameter estimates
#' 2. A vector with labels of parameters
#' 3. Any additional arguments used by the function contained in a single list
#' (called additionalArguments below)
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' CappedL1 regularization:
#' 
#' * Zhang, T. (2010). Analysis of Multi-stage Convex Relaxation for Sparse Regularization. 
#' Journal of Machine Learning Research, 11, 1081–1107.
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param additionalArguments additional argument passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta (see ?controlIsta)
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
#' penModel <- gpCappedL1(
#'   par = b, 
#'   regularized = regularized, 
#'   fn = fittingFunction, 
#'   lambdas = seq(0,100,length.out = 10),
#'   thetas = seq(.1,10,length.out = 10),
#'   additionalArguments = additionalArguments
#' )
#' # optional: plot requires package plotly
#' # plot(penModel)
#' AIC(penModel)
#' @export
gpCappedL1 <- function(par,
                       fn,
                       gr = NULL,
                       ...,
                       regularized,
                       lambdas,
                       thetas,
                       control = controlIsta()){
  
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "cappedL1", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas,
                                                                           alpha = 1), 
                                            method = "ista",
                                            control = control
  )
  
  return(result)
  
}

#' gpLsp 
#' 
#' Implements lsp regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = \lambda \log(1 + |x_j|\theta)}
#' where \eqn{\theta > 0}. 
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take exactly 3 arguments:
#' 
#' 1. A vector with current parameter estimates
#' 2. A vector with labels of parameters
#' 3. Any additional arguments used by the function contained in a single list
#' (called additionalArguments below)
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' lsp regularization:
#' 
#' * Candès, E. J., Wakin, M. B., & Boyd, S. P. (2008). Enhancing Sparsity by 
#' Reweighted l1 Minimization. Journal of Fourier Analysis and Applications, 14(5–6), 
#' 877–905. https://doi.org/10.1007/s00041-008-9045-x
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
gpLsp <- function(par,
                  fn,
                  gr = NULL,
                  ...,
                  regularized,
                  lambdas,
                  thetas,
                  control = controlIsta()){
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "lsp", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas), 
                                            method = "ista",
                                            control = control
  )
  
  return(result)
  
}

#' gpMcp 
#' 
#' Implements mcp regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| - x_j^2/(2\theta) & \text{if } |x_j| \leq \theta\lambda\\
#' \theta\lambda^2/2 & \text{if } |x_j| > \lambda\theta
#' \end{cases}}
#' where \eqn{\theta > 0}. 
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take exactly 3 arguments:
#' 
#' 1. A vector with current parameter estimates
#' 2. A vector with labels of parameters
#' 3. Any additional arguments used by the function contained in a single list
#' (called additionalArguments below)
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' mcp regularization:
#' 
#' * Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax concave penalty. 
#' The Annals of Statistics, 38(2), 894–942. https://doi.org/10.1214/09-AOS729
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
gpMcp <- function(par,
                  fn,
                  gr = NULL,
                  ...,
                  regularized,
                  lambdas,
                  thetas,
                  control = controlIsta()){
  
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "mcp", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas), 
                                            method = "ista",
                                            control = control
  )
  
  return(result)
  
}

#' gpScad 
#' 
#' Implements scad regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| & \text{if } |x_j| \leq \theta\\
#' \frac{-x_j^2 + 2\theta\lambda |x_j| - \lambda^2}{2(\theta -1)} & 
#' \text{if } \lambda < |x_j| \leq \lambda\theta \\
#' (\theta + 1) \lambda^2/2 & \text{if } |x_j| \geq \theta\lambda\\
#' \end{cases}}
#' where \eqn{\theta > 2}.  
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take exactly 3 arguments:
#' 
#' 1. A vector with current parameter estimates
#' 2. A vector with labels of parameters
#' 3. Any additional arguments used by the function contained in a single list
#' (called additionalArguments below)
#' 
#' The gradient function gr is optional. If set to NULL, the \pkg{numDeriv} package
#' will be used to approximate the gradients. Supplying a gradient function
#' can result in considerable speed improvements.  
#' 
#' scad regularization:
#' 
#' * Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized 
#' likelihood and its oracle properties. Journal of the American Statistical Association, 
#' 96(456), 1348–1360. https://doi.org/10.1198/016214501753382273
#'  
#' For more details on GLMNET, see:
#' 
#' * Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
#' Regularization Paths for Generalized Linear Models via Coordinate Descent. 
#' Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
#' * Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
#' A Comparison of Optimization Methods and Software for Large-scale 
#' L1-regularized Linear Classiﬁcation. Journal of Machine Learning Research, 11, 3183–3234.
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
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param additionalArguments additional argument passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' If you are unsure what these parameters are called, use 
#' getLavaanParameters(model) with your lavaan model object
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta (see ?controlIsta)
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
#' penModel <- gpScad(
#'   par = b, 
#'   regularized = regularized, 
#'   fn = fittingFunction, 
#'   lambdas = seq(0,100,length.out = 10),
#'   thetas = seq(2.1,10,length.out = 10),
#'   additionalArguments = additionalArguments
#' )
#' # optional: plot requires package plotly
#' # plot(penModel)
#' AIC(penModel)
#' @export
gpScad <- function(par,
                   fn,
                   gr = NULL,
                   ...,
                   regularized,
                   lambdas,
                   thetas,
                   control = controlIsta()){
  
  removeDotDotDot <- lessSEM::noDotDotDot(fn, ... = ...)
  fn <- removeDotDotDot$fn
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- lessSEM::noDotDotDot(gr, ...)
    gr <- removeDotDotDot$fn
    
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 2)) stop("Theta must be > 2")
  
  result <- lessSEM::gpOptimizationInternal(par = par,
                                            fn = fn,
                                            gr = gr,
                                            additionalArguments = additionalArguments,
                                            penalty = "scad", 
                                            weights = regularized,
                                            tuningParameters = expand.grid(lambda = lambdas, 
                                                                           theta = thetas), 
                                            method = "ista",
                                            control = control
  )
  
  return(result)
}