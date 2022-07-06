test_that("testing general purpose optimization", {
  library(lessSEM)
  library(ncvreg)
  library(glmnet)
  
  ## Testing lasso
  X <- scale(Prostate$X)
  y <- Prostate$y
  
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
  adda <- list(X = X, y = y)
  
  fitFunction <- function(b, parameterLabels, adda){
    N <- length(adda$y)
    pred <- b[1] + adda$X%*%(b[2:length(b)])
    sse <- (.5/N)*sum((adda$y - pred)^2)
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
  testthat::expect_equal(all(round(lassoFit@fits$regM2LL-
                                     penSse,3)==0), TRUE)
  testthat::expect_equal(all(round(lassoFit@fits$m2LL-
                                     sse,3)==0), TRUE)
  
  
  
  mcpFit <- ncvreg(X = X, y = y, penalty = "MCP", nlambda = 30)
  lambdas <- mcpFit$lambda
  thetas <- mcpFit$gamma
  coefs <- coef(mcpFit)
  
  sse <- rep(NA, length(lambdas)) # only one theta was tested
  penSse <- sse
  
  for(l in 1:length(lambdas)){
    N <- length(y)
    sse[l] <- (.5/N)*sum((y - coefs[1,l] - (X %*% coefs[-1,l,drop=FALSE]))^2)
    penSse[l] <- sse[l] + 
      lambdas[l] * sum(abs(coefs[-1,l]))
  }
  
  mcpFitGp <- gpMcp(par = b, 
                    regularized = paste0("b", 1:ncol(X)), 
                    fn = fitFunction, 
                    lambdas = lambdas, 
                    thetas = thetas,
                    additionalArguments = adda, 
                    control = controlIsta())
  
  testthat::expect_equal(all(round(mcpFitGp@fits$regM2LL-
                                     penSse,3)<=0), TRUE)
  
  # plot(mcpFit)
  # plot(mcpFitGp)
  
  scadFit <- ncvreg(X = X, y = y, penalty = "SCAD", nlambda = 30, gamma = 3.4)
  lambdas <- scadFit$lambda
  thetas <- scadFit$gamma
  coefs <- coef(scadFit)
  
  sse <- rep(NA, length(lambdas)) # only one theta was tested
  penSse <- sse
  
  for(l in 1:length(lambdas)){
    N <- length(y)
    sse[l] <- (.5/N)*sum((y - coefs[1,l] - (X %*% coefs[-1,l,drop=FALSE]))^2)
    penSse[l] <- sse[l] + 
      lambdas[l] * sum(abs(coefs[-1,l]))
  }
  
  scadFitGp <- gpScad(par = b, 
                      regularized = paste0("b", 1:ncol(X)), 
                      fn = fitFunction, 
                      lambdas = lambdas, 
                      thetas = thetas,
                      additionalArguments = adda, 
                      control = controlIsta())
  
  testthat::expect_equal(all(round(scadFitGp@fits$regM2LL-
                                     penSse,3)<=0), TRUE)
  
  # plot(scadFit)
  # plot(scadFitGp)
})
