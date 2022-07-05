#' gpCappedL1
#' 
#' This function allows for regularization of custom penatly functions with the
#' cappedL1 penalty.
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
                     additionalArguments = list(),
                     regularized,
                     lambdas,
                     thetas,
                     control = controlIsta()){
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

#' gpScad
#' 
#' This function allows for regularization of custom penatly functions with the
#' scad penalty.
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
                 additionalArguments = list(),
                 regularized,
                 lambdas,
                 thetas,
                 control = controlIsta()){
  
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