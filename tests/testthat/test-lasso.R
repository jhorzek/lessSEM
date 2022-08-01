test_that("testing elasticNet-lasso-c", {
  library(lslx)
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
  
  # fit model using lslx
  lslxModelSyntax <- paste0(paste0('fix(1)*', yNames[1], ' + ', paste0(yNames[2:5], collapse = " + "), " <=: f"),"\n",
                            paste0(paste0(yNames[6:length(yNames)], collapse = " + "), " <~: f"),"\n",
                            paste0(yNames, collapse = " + "), " <= 1"
  )
  fitLslx <- lslx$new(model = lslxModelSyntax,
                      sample_cov = stats::cov(y),
                      sample_size = nrow(y)
  )
  
  #fitLslx$penalize_coefficient(name = paste0("y", 6:ncol(y)," <- f"))
  
  lambdas <- seq(0,.3,.01)
  fitLslx$fit(penalty_method = "lasso",lambda_grid = lambdas, loss = "ml")
  
  # extract fits
  lslxParameter <- matrix(NA, 
                          nrow = length(lambdas), 
                          ncol = length(fitLslx$extract_coefficient(lambda = 0)))
  colnames(lslxParameter) <- names(fitLslx$extract_coefficient(lambda = 0))
  
  m2LLs <- rep(NA, length(lambdas))
  AICs <- rep(NA, length(lambdas))
  BICs <- rep(NA, length(lambdas))
  penalty <- rep(NA, length(lambdas))
  regularized <- paste0("y", 6:ncol(y),"<-f/g") 
  for(l in 1:length(lambdas)){
    pars <- fitLslx$extract_coefficient(lambda = lambdas[l])
    lslxParameter[l,names(pars)] <- pars
    
    # lslx seems to not compute the implied means correctly here:
    #impl_means <- fitLslx$extract_implied_mean(lambda = lambdas[l])$g
    impl_means <- apply(y,2,mean)
    impl_cov <- fitLslx$extract_implied_cov(lambda = lambdas[l])$g
    coefficents <- fitLslx$extract_coefficient(lambda = lambdas[l])
    regularizedCoefficients <- coefficents[regularized]
    m2LLs[l] <- -2*sum(mvtnorm::dmvnorm(x = y, mean = impl_means, sigma = impl_cov, log = TRUE))
    AICs[l] <- m2LLs[l] + 2*(length(coefficents) - sum(regularizedCoefficients == 0) -1) # -1 because lslx also reports the fixed loading
    BICs[l] <- m2LLs[l] + log(N)*(length(coefficents) - sum(regularizedCoefficients == 0) -1)
    penalty[l] <- N*lambdas[l]*sum(abs(regularizedCoefficients))
  }
  
  # replicate with regularizedSEM

  rsemIsta <- lasso(lavaanModel = modelFit, 
                         regularized = paste0("f=~y",6:ncol(y)), 
                         lambdas = lambdas,
                         method = "ista",
                         control = controlIsta(verbose = 0, 
                                               startingValues = "est")
  )
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] - lslxParameter[,regularized]) < .002), TRUE)
  plot(rsemIsta)
  coef(rsemIsta)
  coef(rsemIsta, criterion = "AIC")
  coef(rsemIsta, criterion = "BIC")
  
  testthat::expect_equal(all(round(rsemIsta@fits$regM2LL - (m2LLs + penalty))==0), TRUE)
  testthat::expect_equal(all(round(AIC(rsemIsta)$AIC - (AICs))==0), TRUE)
  testthat::expect_equal(all(round(BIC(rsemIsta)$BIC - (BICs))==0), TRUE)
  
  rsemGlmnet <- lasso(lavaanModel = modelFit, 
                      regularized = paste0("f=~y",6:ncol(y)), 
                           lambdas = lambdas,
                           method = "glmnet",
                           control = controlGlmnet()
  )
  testthat::expect_equal(all(abs(rsemGlmnet@parameters[,rsemGlmnet@regularized] - lslxParameter[,regularized]) < .002), TRUE)
  plot(rsemGlmnet)
  coef(rsemGlmnet)
  coef(rsemGlmnet, criterion = "AIC")
  coef(rsemGlmnet, criterion = "BIC")
  
  testthat::expect_equal(all(round(rsemGlmnet@fits$regM2LL - (m2LLs + penalty))==0), TRUE)
  testthat::expect_equal(all(round(AIC(rsemGlmnet)$AIC - (AICs))==0), TRUE)
  testthat::expect_equal(all(round(BIC(rsemGlmnet)$BIC - (BICs))==0), TRUE)
  
  # set automatic lambda:
  rsem2 <- lessSEM::lasso(lavaanModel = modelFit, 
                                  regularized = paste0("f=~y",6:ncol(y)),
                                  nLambdas = 10)
  testthat::expect_equal(all(apply(rsem2@parameters[,paste0("f=~y",6:ncol(y))] == 0,2,sum) > 0), TRUE)
  
  rsem2 <- lessSEM::lasso(lavaanModel = modelFit, 
                                  regularized = paste0("f=~y",6:ncol(y)),
                                  reverse = FALSE,
                                  nLambdas = 10)
  testthat::expect_equal(all(apply(rsem2@parameters[,paste0("f=~y",6:ncol(y))] == 0,2,sum) > 0), TRUE)
})
