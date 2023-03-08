test_that("testing translation from lessSEM to lavaan", {
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
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  rsem <- lasso(lavaanModel = modelFit, 
                regularized = regularizedLavaan,
                reverse = FALSE,
                nLambdas = 5)
  
  lambdas <- rsem@fits$lambda 
  
  for(l in 1:length(lambdas)){
    lavaanModel <- lessSEM2Lavaan(regularizedSEM = rsem, lambda = lambdas[l])
    testthat::expect_equal(
      all(round(getLavaanParameters(lavaanModel) - unlist(rsem@parameters[l,names(getLavaanParameters(lavaanModel))]),
                10) == 0), 
      TRUE)
  }
  
  rsem <- scad(lavaanModel = modelFit, 
               regularized = regularizedLavaan,
               lambdas = seq(0,1,length.out = 5), 
               thetas = c(2.5,2.6))
  
  lambdas <- rsem@fits$lambda 
  thetas <- rsem@fits$theta 
  
  for(l in 1:length(lambdas)){
      rw <- which(thetas == thetas[l] & lambdas == lambdas[l])
      
      lavaanModel <- lessSEM2Lavaan(regularizedSEM = rsem, lambda = lambdas[l], theta = thetas[l])
      testthat::expect_equal(
        all(round(getLavaanParameters(lavaanModel) - unlist(rsem@parameters[rw,names(getLavaanParameters(lavaanModel))]),
                  10) == 0), 
        TRUE)
  }
  
  testthat::expect_equal(is(lessSEM2Lavaan(regularizedSEM = rsem, criterion = "AIC"), "lavaan"), TRUE)
  
  testthat::expect_equal(is(try(lessSEM2Lavaan(regularizedSEM = rsem, criterion = "GAIC"), silent = TRUE), "try-error"), TRUE)
})
