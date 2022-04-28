test_that("testing approximate influence", {
  
  library(lavaan)
  library(aCV4SEM)
  
  set.seed(123)
  N <- 100
  l1 <- 1; l2 <- .2; l3 <- 0;
  v1 <- .2; v2 <- .8; v3 <- 1
  
  f <- matrix(rnorm(N, 0, 1), ncol = 1)
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
  
  # replicate with regularizedSEM
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  rsem <- regularizeSEM(lavaanModel = modelFit, 
                        regularizedParameterLabels = regularizedLavaan,
                        penalty = "lasso", 
                        nLambdas = 30)
  
  aI <- aI4regularizedSEM(regularizedSEM = rsem, k = N)
  plot(aI)
  AICs <- AIC(aI)
  BICs <- BIC(aI)
  
  coef(AICs)
  coef(BICs)
  
  plot(AICs)
  plot(BICs)
  
  mean(AICs)
  mean(BICs)
  
  plot(mean(AICs)$lambda, mean(AICs)$value,type = "l")
})
