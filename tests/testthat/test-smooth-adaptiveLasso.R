test_that("testing smooth adaptive lasso", {
  # Smooth penalty functions should not be used and are not tested.
  testthat::skip("Skipping test for smooth penalty functions")
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  N <- 100
  l1 <- 1; l2 <- .2; l3 <- 0;
  v1 <- .2; v2 <- .8; v3 <- 1
  
  f <- matrix(stats::rnorm(N, 0, 1), ncol = 1)
  L <- matrix(c(rep(l1,5), rep(l2,5), rep(l3,15)), nrow = 1)
  y <- matrix(NA, nrow = N, ncol = ncol(L))
  
  covs <- c(rep(v1,5), rep(v2,5), rep(v3,15))
  
  for(i in 1:N){
    y[i,] <- L*f[i,] +  mvtnorm::rmvnorm(1, sigma = diag(covs))
  }
  
  yNames <- paste0("y", 1:ncol(y))
  colnames(y) <- yNames
  
  modelSyntax <- paste0('f =~ 1*', yNames[1], ' + ', paste0(yNames[2:length(yNames)], collapse = " + "))
  modelFit = cfa(modelSyntax, y, meanstructure = TRUE)
  
  
  # regularize 
  lambdas <- seq(0,1,.1)
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  rsem <- adaptiveLasso(lavaanModel = modelFit, 
                regularized = regularizedLavaan,
                lambdas = lambdas)
  
  
  
  apprRegsem <- smoothAdaptiveLasso(lavaanModel = modelFit, 
                          regularized = regularizedLavaan,
                          epsilon = 1e-8,
                          tau = 1e-4,
                          lambdas = lambdas)
  
  parameterDifference <- apprRegsem@parameters[,rsem@parameterLabels] - rsem@parameters[,rsem@parameterLabels]
  matplot(abs(parameterDifference[,regularizedLavaan]), type  ="l")
  testthat::expect_equal(max(abs(parameterDifference[,regularizedLavaan])) < .03, TRUE)
})
