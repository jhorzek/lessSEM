test_that("testing maxLambda", {
  library(regsem)
  library(lessSEM)
  set.seed(123)
  N <- 50
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
  
  #### Regularize ####
  lavaanParameters <- lessSEM::getLavaanParameters(modelFit)
  weights <- rep(0, length(lavaanParameters))
  names(weights) <- names(lavaanParameters)
  weights[paste0("f=~y",6:ncol(y))] <- 1/abs(lavaanParameters[paste0("f=~y",6:ncol(y))])
  
  rsemIsta <- lasso(lavaanModel = modelFit, 
                         regularized = paste0("f=~y",6:ncol(y)), 
                         nLambdas = 30,
                         method = "ista",
                         control = controlIsta()
  )
  testthat::expect_equal(any(apply(rsemIsta@parameters[,rsemIsta@regularized] == 0, 1, all)), TRUE)
  
  plot(rsemIsta)
  coef(rsemIsta)
  coef(rsemIsta, criterion = "AIC")
  coef(rsemIsta, criterion = "BIC")
  
  rsemGlmnet <- lasso(lavaanModel = modelFit, 
                           regularized = paste0("f=~y",6:ncol(y)), 
                           nLambdas = 30,
                           method = "glmnet",
                           control = controlGlmnet()
  )
  testthat::expect_equal(any(apply(rsemGlmnet@parameters[,rsemGlmnet@regularized] == 0, 1, all)), TRUE)
  
  plot(rsemGlmnet)
  coef(rsemGlmnet)
  coef(rsemGlmnet, criterion = "AIC")
  coef(rsemGlmnet, criterion = "BIC")
  
  
})
