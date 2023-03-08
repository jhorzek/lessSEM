test_that("testing new tau", {
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
                nLambdas = 5)
  
  lambdas <- rsem@fits$lambda 
  
  apprRegsem <- smoothLasso(lavaanModel = modelFit, 
                            regularized = regularizedLavaan, 
                            epsilon = 1e-8, 
                            tau = 0,
                            lambdas = lambdas)
  
  testthat::expect_equal(length(unique(apprRegsem@fits$nonZeroParameters)) , 1)
  
  newTaus <- seq(1e-10,1e-1,length.out = 10)
  for(nt in newTaus){
    apprRegsem <- newTau(apprRegsem, nt)
    nZero <- apply(apprRegsem@parameters[,apprRegsem@regularized], 1, function(x) sum(abs(x) <= nt))
    nZero[which(apprRegsem@fits$lambda == 0)] <- 0
    testthat::expect_equal(all(apprRegsem@fits$nonZeroParameters == (length(apprRegsem@parameterLabels) - nZero)) , TRUE)
  }
})
