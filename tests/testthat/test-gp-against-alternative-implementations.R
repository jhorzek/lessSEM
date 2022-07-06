test_that("testing general purpose optimization", {
  library(lessSEM)
  library(ncvreg)
  library(glmnet)
  
  ## Testing lasso
  X <- scale(Prostate$X)
  y <- Prostate$y
  N <- length(y)
  
  out <- glmnet(x = X, 
                y = y, 
                lambda = seq(0,.8,.1),
                standardize = FALSE)
  lambdas <- out$lambda
  
  coefs <- coef(out)
  sse <- rep(NA, length(lambdas))
  penSse <- sse
  for(l in 1:length(lambdas)){
    sse[l] <- sum((y - coefs[1,l] - (X %*% coefs[-1,l,drop=FALSE]))^2)
    penSse[l] <- sse[l] + 
      lambdas[l] * sum(abs(coefs[-1,l]))
  }
  obj_function(y, coefs[1,l] + (X %*% coefs[-1,l,drop=FALSE]), weights, family, lambda, alpha, coefficients, vp)
  
  
  
  # with lessSEM
  adda <- list(X = X, y = y)
  
  fitFunction <- function(b, parameterLabels, adda){
    pred <- b[1] + adda$X%*%(b[2:length(b)])
    sse <- sum((adda$y - pred)^2)
    return(sse)
  }
  
  # initialize
  b <- rep(0, ncol(X)+1)
  names(b) <- paste0("b",0:(length(b)-1))
  fitFunction(matrix(b,nrow = 1), names(b), adda)
  
  lassoFit <- gpLasso(par = b, 
                      regularized = paste0("b", 1:ncol(X)), 
                      fn = fitFunction, 
                      lambdas = lambdas, 
                      additionalArguments = adda, 
                      method = "glmnet", 
                      control = controlGlmnet())
  penSse
  lassoFit@fits$regM2LL
  sse 
  lassoFit@fits$m2LL
  
  l <- 3
  bs <- fit$beta[,l,drop = FALSE]
  sum((y - bs[1] - (X %*% bs[-1,]))^2) + 
    lambdas[l] * sum(abs(bs[-1,]))
  fits[l]
  
  bs <- coefs[,l,drop = FALSE]
  sum((y - bs[1] - (X %*% bs[-1,]))^2) + 
    lambdas[l] * sum(abs(bs[-1,]))
  fits[l]
  
  
  bs2 <- unlist(lassoFit@parameters[l,lassoFit@parameterLabels])
  sum((y - cbind(1,fit$X)%*%bs2)^2) + 
    lambdas[l] * sum(abs(bs2[2:length(bs2)]))
  lassoFit@fits$regM2LL[l]
  
  
  lassoFit@fits$m2LL
  
  plot(fit)
  plot(lassoFit)
  
  fitFunctionLasso <- function(b, parameterLabels, adda){
    pred <- b[1] + adda$X%*%(b[2:length(b)])
    sse <- sum((adda$y - pred)^2)
    penFit <- sse + adda$lambda*sum(sqrt(b[2:length(b)]^2 + 1e-10))
    return(penFit)
  }
  
  for(l in 1:length(lambdas)){
    adda$lambda <- lambdas[l]
    
    print("Optim:")
    opt <- optim(par = b, 
                 fn = fitFunctionLasso, 
                 parameterLabels = names(b), 
                 adda = adda, 
                 method = "BFGS")
    print(
      opt$value
    )
    print(fitFunctionLasso(b = opt$par, parameterLabels = names(b), adda = adda))
    
    print("gp:")
    print(
      lassoFit@fits$regM2LL[l]
    )
    
    print("lasso:")
    bs <- fit$beta[,l]
    print(fitFunctionLasso(b = bs, parameterLabels = names(b), adda = adda))
    print(fit$loss[l])
  }
  
  # initialize
  b <- rep(0, ncol(X)+1)
  names(b) <- paste0("b",0:(length(b)-1))
  fitFunction(matrix(b,nrow = 1), names(b), adda)
  
  optim(par = b, fn = fitFunction, parameterLabels = names(b), )
  
  
  
  for(lambda in lambdas)
    testthat::expect_equal(all(round(coef(fit, lambda = lambdas[1]) -
                                       coef(lassoFit)[coef(lassoFit)$lambda == length(y)*lambdas[1],names(b)],
                                     4) == 0),
                           TRUE)
  
  control <- list(
    initialHessian = diag(1,length(b)),
    stepSize = .9,
    sigma = 0,
    gamma = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 1000,
    breakOuter = 1e-6,
    breakInner = 1e-5,
    convergenceCriterion = 0, 
    verbose = 0
  )
  IL <- new(glmnetEnetGeneralPurpose, 
            weights, 
            control)
  
  lassoResult <- IL$optimize(
    b,
    fitFunction,
    gradientFunction,
    listElements,
    500,
    1
  )
  lassoResult$fit
  lassoResult$rawParameters
  coef(lmFit)
  testthat::expect_equal(all(round(lassoResult$rawParameters[c("b4", "b5")] -
                                     coef(lm(y ~ 0+x4+x5, df)),
                                   4) == 0),
                         TRUE)
  
})
