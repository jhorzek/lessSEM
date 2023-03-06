test_that("testing general purpose optimization with gr", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmnet")
  testthat::skip_if_not_installed("ncvreg")
  library(lessSEM)
  library(glmnet)
  library(ncvreg)
  set.seed(123)
  
  # first, we simulate data for our
  # linear regression.
  N <- 100 # number of persons
  p <- 10 # number of predictors
  X <- matrix(stats::rnorm(N*p),	nrow = N, ncol = p) # design matrix
  b <- c(rep(1,4),
         rep(0,6)) # true regression weights
  y <- X%*%matrix(b,ncol = 1) + stats::rnorm(N,0,.2)
  
  # First, we must construct a fiting function
  # which returns a single value. We will use
  # the residual sum squared as fitting function.
  
  # Let's start setting up the fitting function:
  sseFun <- function(par, y, X, N){
    # par is the parameter vector
    # y is the observed dependent variable
    # X is the design matrix
    # N is the sample size
    pred <- X %*% matrix(par, ncol = 1) #be explicit here:
    # we need par to be a column vector
    sse <- sum((y - pred)^2)
    # we scale with .5/N to get the same results as glmnet
    return((.5/N)*sse)
  }
  
  sseGrad <- function(par, y, X, N){
    
      gradients = (-2.0*t(X) %*% y + 2.0*t(X)%*%X%*%matrix(par,ncol = 1))
      
      gradients = (.5/length(y))*gradients
      return(t(gradients))
  }
  
  # add intercept
  Xext <- cbind(1,X)
  
  # let's define the starting values:
  b <- c(solve(t(Xext)%*%Xext)%*%t(Xext)%*%y) # we will use the lm estimates
  names(b) <- paste0("b", 1:length(b))
  # names of regularized parameters
  regularized <- paste0("b",2:length(b))
  
  sseGrad(b,y,Xext,N)
  
  # optimize
  ridgePen <- gpRidge(
    par = b,
    regularized = regularized,
    fn = sseFun,
    gr = sseGrad,
    lambdas = seq(0,1,.1),
    X = Xext,
    y = y,
    N = N
  )
  
  # for comparison:
  fittingFunction <- function(par, y, X, N, lambda){
    pred <- X %*% matrix(par, ncol = 1)
    sse <- sum((y - pred)^2)
    return((.5/N)*sse + lambda * sum(par[2:length(par)]^2))
  }
  
  testthat::expect_equal(
    all(round(optim(par = b,
                    fn = fittingFunction,
                    y = y,
                    X = Xext,
                    N = N,
                    lambda =  ridgePen@fits$lambda[3],
                    method = "BFGS")$par - 
                ridgePen@parameters[3,ridgePen@parameterLabels], 4) == 0), TRUE)
  
  ## Testing lasso
  X <- scale(X)
  out <- glmnet(x = X, 
                y = y, 
                lambda = seq(0,.8,.1),
                standardize = FALSE)
  lambdas <- out$lambda
  
  coefs <- coef(out)
  sse <- rep(NA, length(lambdas))
  penSse <- sse
  for(l in 1:length(lambdas)){
    N <- length(y)
    sse[l] <- (.5/N)*sum((y - coefs[1,l] - (X %*% coefs[-1,l,drop=FALSE]))^2)
    penSse[l] <- sse[l] + 
      lambdas[l] * sum(abs(coefs[-1,l]))
  }
  
  # with lessSEM
  lassoPen <- gpLasso(
    par = b, 
    regularized = regularized, 
    fn = sseFun, 
    gr = sseGrad,
    lambdas = lambdas, 
    X = cbind(1,X),
    y = y,
    N = N
  )
  
  testthat::expect_equal(all(abs(t(lassoPen@parameters[,lassoPen@parameterLabels]) -
                                   coef(out)) <.0001), TRUE)
  
  testthat::expect_equal(all(abs(lassoPen@fits$regM2LL-
                                   penSse) <.0001), TRUE)
  
  
  X <- ncvreg::std(X)
  attr(X, "center") <- NULL
  attr(X, "scale") <- NULL
  attr(X, "nonsingular") <- NULL
  
  mcpFit <- ncvreg(X = X, y = y, penalty = "MCP", nlambda = 30)
  lambdas <- mcpFit$lambda
  thetas <- mcpFit$gamma
  coefs <- coef(mcpFit)
  
  mcpFitGp <- gpMcp(par = b, 
                    regularized = regularized, 
                    fn = sseFun, 
                    gr = sseGrad,
                    lambdas = lambdas, 
                    thetas = thetas,
                    X = cbind(1,X),
                    y = y,
                    N = N)
  
  testthat::expect_equal(all(abs(t(mcpFitGp@parameters[,mcpFitGp@parameterLabels]) -
                                   coef(mcpFit)) <.001), TRUE)
  
  # plot(mcpFit)
  # plot(mcpFitGp)
  
  scadFit <- ncvreg(X = X, y = y, penalty = "SCAD", nlambda = 30, gamma = 3.4)
  lambdas <- scadFit$lambda
  thetas <- scadFit$gamma
  coefs <- coef(scadFit)
  
  scadFitGp <- gpScad(par = b, 
                      regularized = regularized, 
                      fn = sseFun, 
                      gr = sseGrad,
                      lambdas = lambdas, 
                      thetas = thetas,
                      X = cbind(1,X),
                      y = y,
                      N = N)
  
  testthat::expect_equal(all(abs(t(scadFitGp@parameters[,scadFitGp@parameterLabels]) -
                                   coef(scadFit)) <.001), TRUE)
  
  # plot(scadFit)
  # plot(scadFitGp)
})

