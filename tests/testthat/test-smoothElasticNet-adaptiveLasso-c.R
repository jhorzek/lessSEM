test_that("testing smoothElasticNet-ridge-c", {
  library(lslx)
  library(lavaan)
  library(linr)
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
  
  # fit model using lslx
  lslxModelSyntax <- paste0('fix(1)*', yNames[1], ' + ', paste0(yNames[2:length(yNames)], collapse = " + "), " <=: f") 
  fitLslx <- lslx$new(model = lslxModelSyntax,
                      sample_cov = cov(y),
                      sample_size = nrow(y)
  )
  
  fitLslx$penalize_coefficient(name = paste0("y", 6:ncol(y)," <- f"))
  
  lambdas <- seq(0,.8,.05)
  fitLslx$fit(penalty_method = "ridge",lambda_grid = lambdas, loss = "ml")
  
  # extract fits
  lslxParameter <- matrix(NA, 
                          nrow = length(lambdas), 
                          ncol = length(fitLslx$extract_coefficient(lambda = 0)))
  colnames(lslxParameter) <- names(fitLslx$extract_coefficient(lambda = 0))
  
  for(l in 1:length(lambdas)){
    pars <- fitLslx$extract_coefficient(lambda = lambdas[l])
    lslxParameter[l,names(pars)] <- pars
  }
  regularized <- paste0("y", 6:ncol(y),"<-f/g") 
  
  # replicate with regularizedSEM
  lavaanParameters <- linr::getLavaanParameters(modelFit)
  weights <- rep(0, length(lavaanParameters))
  names(weights) <- names(lavaanParameters)
  weights[paste0("f=~y",6:ncol(y))] <- 1
  
  rsemBfgs <- smoothElasticNetValue(lavaanModel = modelFit, 
                               weights =  weights, 
                               alphas = 0, 
                               lambdas = lambdas,
                               epsilon = 1e-8, 
                               control = controlBFGS(startingValues = "start", 
                                                     initialHessian = diag(100, length(lavaanParameters)),
                                                     breakOuter = 1e-5, 
                                                     breakInner = 1e-5, 
                                                     sigma = 0,
                                                     verbose = 0)
  )
  rsemBfgs@fits$regM2LL
  
  testthat::expect_equal(all(abs(rsemBfgs@parameters[,rsemBfgs@regularized] - lslxParameter[,regularized]) < .02), TRUE)
  plot(rsemBfgs)
  coef(rsemBfgs)
  coef(rsemBfgs, alpha = 0, lambda = .1)
  
  rsemGlmnet <- elasticNet(lavaanModel = modelFit, 
                           weights =  weights, 
                           alphas = 0, 
                           lambdas = lambdas,
                           method = "glmnet",
                           control = controlGlmnet(verbose = 0, 
                                                   startingValues = "est")
  )
  testthat::expect_equal(all(abs(rsemGlmnet@parameters[,rsemGlmnet@regularized] - rsemBfgs@parameters[,rsemBfgs@regularized]) < .02), TRUE)
  plot(rsemGlmnet)
  coef(rsemGlmnet)
  coef(rsemGlmnet, alpha = 0, lambda = .1)
  
  ## Test exact cross-validation
  warning("Not testing approximate cross-validation")
  # cvExact <- CV4regularizedSEM(regularizedSEM = rsem, k = N)
  # coef(cvExact)
  # coef(cvExact, rule = "1sd")
  # coef(cvExact, rule = "penalized")
  # coef(cvExact, alpha = 1, lambda = .1)
  # plot(cvExact)
  
})
