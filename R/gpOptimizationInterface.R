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
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param reverse if set to TRUE and nLambdas is used, lessSEM will start with the
#' largest lambda and gradually decrease lambda. Otherwise, lessSEM will start with
#' the smallest lambda and gradually increase it.
#' @param curve Allows for unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). 
#' If curve is close to 1 all lambda values will be equally spaced, if curve is large 
#' lambda values will be more concentrated close to 0. See ?lessSEM::curveLambda for more information.
#' @param ... additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

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
#' 
#' # You can access the fit results as follows:
#' lassoPen@fits
#' # Note that we won't compute any fit measures automatically, as
#' # we cannot be sure how the AIC, BIC, etc are defined for your objective function 
#' 
#' @export
gpLasso <- function(par,
                    regularized,
                    fn,
                    gr = NULL,
                    lambdas = NULL,
                    nLambdas = NULL,
                    reverse = TRUE,
                    curve = 1,
                    ...,
                    method = "glmnet", 
                    control = lessSEM::controlGlmnet()
){
  
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(is.null(lambdas) && is.null(nLambdas)){
    stop("Specify either lambdas or nLambdas")
  }
  
  if(!is.null(nLambdas)){
    tuningParameters <- data.frame(nLambdas = nLambdas,
                                   reverse = reverse,
                                   curve = curve)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1,
                                   theta = 0)
  }
  
  result <- .gpOptimizationInternal(par = par,
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
#' function. This fitting functions _must_ take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param weights labeled vector with adaptive lasso weights. NULL will use 1/abs(par)
#' @param fn R function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr R function which takes the parameters as input and returns the 
#' gradients of the objective function. If set to NULL, numDeriv will be used
#' to approximate the gradients 
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
#' @param ... additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

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
#' b <- c(solve(t(X)%*%X)%*%t(X)%*%y) # we will use the lm estimates
#' names(b) <- paste0("b", 1:length(b))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # define the weight for each of the parameters
#' weights <- 1/abs(b)
#' # we will re-scale the weights for equivalence to glmnet.
#' # see ?glmnet for more details
#' weights <- length(b)*weights/sum(weights)
#' 
#' # optimize
#' adaptiveLassoPen <- gpAdaptiveLasso(
#'   par = b, 
#'   regularized = regularized, 
#'   weights = weights,
#'   fn = fittingFunction, 
#'   lambdas = seq(0,1,.01), 
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' plot(adaptiveLassoPen)
#' # You can access the fit results as follows:
#' adaptiveLassoPen@fits
#' # Note that we won't compute any fit measures automatically, as
#' # we cannot be sure how the AIC, BIC, etc are defined for your objective function 
#' 
#' # for comparison:
#' # library(glmnet)
#' # coef(glmnet(x = X,
#' #            y = y,
#' #            penalty.factor = weights,
#' #            lambda = adaptiveLassoPen@fits$lambda[20],
#' #            intercept = FALSE,
#' #            standardize = FALSE))[,1]
#' # adaptiveLassoPen@parameters[20,]
#' @export
gpAdaptiveLasso <- function(par,
                            regularized,
                            weights = NULL,
                            fn,
                            gr = NULL,
                            lambdas = NULL,
                            nLambdas = NULL,
                            reverse = TRUE,
                            curve = 1,
                            ...,
                            method = "glmnet", 
                            control = lessSEM::controlGlmnet()){
  
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(is.null(weights)){
    weights <- 1/abs(par)
    weights[!names(weights) %in% regularized] <- 0
    cat("\n")
    rlang::inform(c("Note","Building weights based on par as weights = 1/abs(par)."))
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
    tuningParameters <- data.frame(nLambdas = nLambdas,
                                   reverse = reverse,
                                   curve = curve)
  }else{
    tuningParameters <- data.frame(lambda = lambdas,
                                   alpha = 1,
                                   theta = 0)
  }
  
  result <- .gpOptimizationInternal(par = par,
                                    fn = fn,
                                    gr = gr,
                                    additionalArguments = additionalArguments,
                                    penalty = "adaptiveLasso", 
                                    weights = weights,
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
#' function. This fitting functions _must_ take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters as input and returns the 
#' fit value (a single value)
#' @param gr R function which takes the parameters as input and returns the 
#' gradients of the objective function. If set to NULL, numDeriv will be used
#' to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param ... additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

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
#' b <- c(solve(t(X)%*%X)%*%t(X)%*%y) # we will use the lm estimates
#' names(b) <- paste0("b", 1:length(b))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # optimize
#' ridgePen <- gpRidge(
#'   par = b,
#'   regularized = regularized,
#'   fn = fittingFunction,
#'   lambdas = seq(0,1,.01),
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' plot(ridgePen)
#' 
#' # for comparison:
#' # fittingFunction <- function(par, y, X, N, lambda){
#' #   pred <- X %*% matrix(par, ncol = 1) 
#' #   sse <- sum((y - pred)^2)
#' #   return((.5/N)*sse + lambda * sum(par^2))
#' # }
#' # 
#' # optim(par = b, 
#' #       fn = fittingFunction, 
#' #       y = y,
#' #       X = X,
#' #       N = N,
#' #       lambda =  ridgePen@fits$lambda[20], 
#' #       method = "BFGS")$par
#' # ridgePen@parameters[20,]
#' @export
gpRidge <- function(par,
                    regularized,
                    fn,
                    gr = NULL,
                    lambdas,
                    ...,
                    method = "glmnet", 
                    control = lessSEM::controlGlmnet()){
  
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  
  tuningParameters <- data.frame(lambda = lambdas,
                                 alpha = 0,
                                 theta = 0)
  
  result <- .gpOptimizationInternal(par = par,
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
#' function. This fitting functions _must_ take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param par labeled vector with starting values
#' @param regularized vector with names of parameters which are to be regularized.
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param alphas numeric vector with values of the tuning parameter alpha. Must be
#' between 0 and 1. 0 = ridge, 1 = lasso.
#' @param regularized vector with names of parameters which are to be regularized.
#' @param ... additional arguments passed to fn and gr
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized
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
#' b <- c(solve(t(X)%*%X)%*%t(X)%*%y) # we will use the lm estimates
#' names(b) <- paste0("b", 1:length(b))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # optimize
#' elasticNetPen <- gpElasticNet(
#'   par = b,
#'   regularized = regularized,
#'   fn = fittingFunction,
#'   lambdas = seq(0,1,.1),
#'   alphas = c(0, .5, 1),
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' 
#' # optional: plot requires plotly package
#' # plot(elasticNetPen)
#' 
#' # for comparison:
#' fittingFunction <- function(par, y, X, N, lambda, alpha){
#'   pred <- X %*% matrix(par, ncol = 1)
#'   sse <- sum((y - pred)^2)
#'   return((.5/N)*sse + (1-alpha)*lambda * sum(par^2) + alpha*lambda *sum(sqrt(par^2 + 1e-8)))
#' }
#' 
#' round(
#'   optim(par = b,
#'       fn = fittingFunction,
#'       y = y,
#'       X = X,
#'       N = N,
#'       lambda =  elasticNetPen@fits$lambda[15],
#'       alpha =  elasticNetPen@fits$alpha[15],
#'       method = "BFGS")$par,
#'   4)
#' elasticNetPen@parameters[15,]
#' @export
gpElasticNet <- function(par,
                         regularized,
                         fn,
                         gr = NULL,
                         lambdas,
                         alphas,
                         ...,
                         method = "glmnet", 
                         control = lessSEM::controlGlmnet()){
  
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  tuningParameters <- expand.grid(lambda = lambdas,
                                  alpha = alphas,
                                  theta = 0)
  
  
  result <- .gpOptimizationInternal(par = par,
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
#' function. This fitting functions _must_ take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param ... additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas parameters whose absolute value is above this threshold will be penalized with
#' a constant (theta)
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

#' @examples 
#' # This example shows how to use the optimizers
#' # for other objective functions. We will use
#' # a linear regression as an example. Note that
#' # this is not a useful application of the optimizers
#' # as there are specialized packages for linear regression
#' # (e.g., glmnet)
#' 
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
#' b <- c(solve(t(X)%*%X)%*%t(X)%*%y) # we will use the lm estimates
#' names(b) <- paste0("b", 1:length(b))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # optimize
#' cL1 <- gpCappedL1(
#'   par = b,
#'   regularized = regularized,
#'   fn = fittingFunction,
#'   lambdas = seq(0,1,.1),
#'   thetas = c(0.001, .5, 1),
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' 
#' # optional: plot requires plotly package
#' # plot(cL1)
#' 
#' # for comparison
#' 
#' fittingFunction <- function(par, y, X, N, lambda, theta){
#'   pred <- X %*% matrix(par, ncol = 1)
#'   sse <- sum((y - pred)^2)
#'   smoothAbs <- sqrt(par^2 + 1e-8)
#'   pen <- lambda * ifelse(smoothAbs < theta, smoothAbs, theta)
#'   return((.5/N)*sse + sum(pen))
#' }
#' 
#' round(
#'   optim(par = b,
#'       fn = fittingFunction,
#'       y = y,
#'       X = X,
#'       N = N,
#'       lambda =  cL1@fits$lambda[15],
#'       theta =  cL1@fits$theta[15],
#'       method = "BFGS")$par,
#'   4)
#' cL1@parameters[15,]
#' @export
gpCappedL1 <- function(par,
                       fn,
                       gr = NULL,
                       ...,
                       regularized,
                       lambdas,
                       thetas,
                       method = "glmnet", 
                       control = lessSEM::controlGlmnet()){
  
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .gpOptimizationInternal(par = par,
                                    fn = fn,
                                    gr = gr,
                                    additionalArguments = additionalArguments,
                                    penalty = "cappedL1", 
                                    weights = regularized,
                                    tuningParameters = expand.grid(lambda = lambdas, 
                                                                   theta = thetas,
                                                                   alpha = 1), 
                                    method = method,
                                    control = control
  )
  
  return(result)
  
}

#' gpLsp 
#' 
#' Implements lsp regularization for general purpose optimization problems.
#' The penalty function is given by:
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector must have labels) and a fitting
#' function. This fitting functions must take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param ... additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas numeric vector: values for the tuning parameter theta
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized
#' @examples 
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
#' b <- c(solve(t(X)%*%X)%*%t(X)%*%y) # we will use the lm estimates
#' names(b) <- paste0("b", 1:length(b))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # optimize
#' lspPen <- gpLsp(
#'   par = b,
#'   regularized = regularized,
#'   fn = fittingFunction,
#'   lambdas = seq(0,1,.1),
#'   thetas = c(0.001, .5, 1),
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' 
#' # optional: plot requires plotly package
#' # plot(lspPen)
#' 
#' # for comparison
#' 
#' fittingFunction <- function(par, y, X, N, lambda, theta){
#'   pred <- X %*% matrix(par, ncol = 1)
#'   sse <- sum((y - pred)^2)
#'   smoothAbs <- sqrt(par^2 + 1e-8)
#'   pen <- lambda * log(1.0 + smoothAbs / theta)
#'   return((.5/N)*sse + sum(pen))
#' }
#' 
#' round(
#'   optim(par = b,
#'       fn = fittingFunction,
#'       y = y,
#'       X = X,
#'       N = N,
#'       lambda =  lspPen@fits$lambda[15],
#'       theta =  lspPen@fits$theta[15],
#'       method = "BFGS")$par,
#'   4)
#' lspPen@parameters[15,]
#' @export
gpLsp <- function(par,
                  fn,
                  gr = NULL,
                  ...,
                  regularized,
                  lambdas,
                  thetas,
                  method = "glmnet", 
                  control = lessSEM::controlGlmnet()){
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .gpOptimizationInternal(par = par,
                                    fn = fn,
                                    gr = gr,
                                    additionalArguments = additionalArguments,
                                    penalty = "lsp", 
                                    weights = regularized,
                                    tuningParameters = expand.grid(lambda = lambdas, 
                                                                   theta = thetas,
                                                                   alpha = 1), 
                                    method = method,
                                    control = control
  )
  
  return(result)
  
}

#' gpMcp 
#' 
#' Implements mcp regularization for general purpose optimization problems.
#' The penalty function is given by:
#' \ifelse{html}{\deqn{p( x_j) = \begin{cases}
#' \lambda |x_j| - x_j^2/(2\theta) & \text{if } |x_j| \leq \theta\lambda\\
#' \theta\lambda^2/2 & \text{if } |x_j| > \lambda\theta
#' \end{cases}} where \eqn{\theta > 0}.}{
#' Equation Omitted in Pdf Documentation.}
#' 
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param ... additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas numeric vector: values for the tuning parameter theta
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

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
#' # first, let's add an intercept
#' X <- cbind(1, X)
#' 
#' b <- c(solve(t(X)%*%X)%*%t(X)%*%y) # we will use the lm estimates
#' names(b) <- paste0("b", 0:(length(b)-1))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # optimize
#' mcpPen <- gpMcp(
#'   par = b,
#'   regularized = regularized,
#'   fn = fittingFunction,
#'   lambdas = seq(0,1,.1),
#'   thetas = c(1.001, 1.5, 2),
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' 
#' # optional: plot requires plotly package
#' # plot(mcpPen)
#' 
#' @export
gpMcp <- function(par,
                  fn,
                  gr = NULL,
                  ...,
                  regularized,
                  lambdas,
                  thetas,
                  method = "glmnet", 
                  control = lessSEM::controlGlmnet()){
  
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 0)) stop("Theta must be > 0")
  
  result <- .gpOptimizationInternal(par = par,
                                    fn = fn,
                                    gr = gr,
                                    additionalArguments = additionalArguments,
                                    penalty = "mcp", 
                                    weights = regularized,
                                    tuningParameters = expand.grid(lambda = lambdas, 
                                                                   theta = thetas,
                                                                   alpha = 1), 
                                    method = method,
                                    control = control
  )
  
  return(result)
  
}

#' gpScad 
#' 
#' Implements scad regularization for general purpose optimization problems.
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
#' The interface is similar to that of optim. Users have to supply a vector 
#' with starting values (important: This vector _must_ have labels) and a fitting
#' function. This fitting functions _must_ take a labeled vector with parameter
#' values as first argument. The remaining arguments are passed with the ... argument.
#' This is similar to optim.
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
#' @param par labeled vector with starting values
#' @param fn R function which takes the parameters AND their labels 
#' as input and returns the fit value (a single value)
#' @param gr R function which takes the parameters AND their labels
#' as input and returns the gradients of the objective function. 
#' If set to NULL, numDeriv will be used to approximate the gradients 
#' @param ... additional arguments passed to fn and gr
#' @param regularized vector with names of parameters which are to be regularized.
#' @param lambdas numeric vector: values for the tuning parameter lambda
#' @param thetas numeric vector: values for the tuning parameter theta
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. 
#' @param control used to control the optimizer. This element is generated with 
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet
#' for more details.
#' @returns Object of class gpRegularized

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
#' # first, let's add an intercept
#' X <- cbind(1, X)
#' 
#' b <- c(solve(t(X)%*%X)%*%t(X)%*%y) # we will use the lm estimates
#' names(b) <- paste0("b", 0:(length(b)-1))
#' # names of regularized parameters
#' regularized <- paste0("b",1:p)
#' 
#' # optimize
#' scadPen <- gpScad(
#'   par = b,
#'   regularized = regularized,
#'   fn = fittingFunction,
#'   lambdas = seq(0,1,.1),
#'   thetas = c(2.001, 2.5, 5),
#'   X = X,
#'   y = y,
#'   N = N
#' )
#' 
#' # optional: plot requires plotly package
#' # plot(scadPen)
#' 
#' # for comparison
#' #library(ncvreg)
#' #scadFit <- ncvreg(X = X[,-1], 
#' #                  y = y, 
#' #                  penalty = "SCAD",
#' #                  lambda =  scadPen@fits$lambda[15],
#' #                  gamma =  scadPen@fits$theta[15])
#' #coef(scadFit)
#' #scadPen@parameters[15,]
#' @export
gpScad <- function(par,
                   fn,
                   gr = NULL,
                   ...,
                   regularized,
                   lambdas,
                   thetas,
                   method = "glmnet", 
                   control = lessSEM::controlGlmnet()){
  
  removeDotDotDot <- .noDotDotDot(fn, fnName = "fn", ... = ...)
  fn <- removeDotDotDot[[1]]
  additionalArguments <- removeDotDotDot$additionalArguments
  
  if(!is.null(gr)){
    
    removeDotDotDot <- .noDotDotDot(gr, fnName = "gr", ...)
    gr <- removeDotDotDot[[1]]
    additionalArguments <- c(additionalArguments, 
                             removeDotDotDot$additionalArguments[!names(removeDotDotDot$additionalArguments) %in% names(additionalArguments)])
  }
  
  # remove the ... stuff so that it does not interfere with anything else
  rm("...")
  
  if(any(thetas <= 2)) stop("Theta must be > 2")
  
  result <- .gpOptimizationInternal(par = par,
                                    fn = fn,
                                    gr = gr,
                                    additionalArguments = additionalArguments,
                                    penalty = "scad", 
                                    weights = regularized,
                                    tuningParameters = expand.grid(lambda = lambdas, 
                                                                   theta = thetas,
                                                                   alpha = 1), 
                                    method = method,
                                    control = control
  )
  
  return(result)
}