test_that("testing adaptive lasso", {
  library(regsem)
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
  
  # fit model using regsem
  pars_pen = 5:24
  regsem_cvFit <- cv_regsem(model = modelFit, 
                            pars_pen = pars_pen, 
                            type = "alasso", 
                            n.lambda = 5,
                            gradFun = "ram")
  plot(regsem_cvFit)
  
  regsemPars <- linr::cvregsem2LavaanParameters(cvregsemModel = regsem_cvFit, 
                                                lavaanModel = modelFit)
  
  ## reconstruct fit for regsem
  N <- nrow(y)
  regsemNonConverged <- (regsem_cvFit$fits[,"conv"] != 0) | (regsem_cvFit$fits[,"chisq"] < 0)
  
  saturatedM2LL <- -2*sum(mvtnorm::dmvnorm(
    x = y, 
    mean = apply(y,2,mean), 
    sigma = ((N - 1)/N)*cov(y),
    log = TRUE
  ))
  regsemChisq <- regsem_cvFit$fits[,"chisq"]
  regsemM2LL <- regsemChisq + saturatedM2LL
  
  # replicate with regularizedSEM
  lambdas <- regsem_cvFit$fits[,"lambda"]
  lavaanParameters <- linr::getLavaanParameters(modelFit)
  weights <- rep(0, length(lavaanParameters))
  names(weights) <- names(lavaanParameters)
  weights[paste0("f=~y",6:ncol(y))] <- 1/abs(lavaanParameters[paste0("f=~y",6:ncol(y))])
  
  rsemIsta <- elasticNet(lavaanModel = modelFit, 
                         weights =  weights, 
                         alphas = 1, 
                         lambdas = lambdas,
                         method = "ista",
                         control = controlIsta(verbose = 0, 
                                               startingValues = "est")
  )
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,names(lavaanParameters)] -
                                   regsemPars[,names(lavaanParameters)]) < .1), TRUE)
  plot(rsemIsta)
  coef(rsemIsta)
  coef(rsemIsta, alpha = 1, lambda = lambdas[1])
  coef(rsemIsta, criterion = "AIC")
  coef(rsemIsta, criterion = "BIC")
  
  testthat::expect_equal(all(round(rsemIsta@fits$m2LL - (regsemM2LL))==0), TRUE)
  penalty <- N*lambdas*apply(regsemPars[,names(lavaanParameters)],1, function(x) sum(abs(x)*weights))
  testthat::expect_equal(all(round(rsemIsta@fits$regM2LL - (regsemM2LL+penalty))==0), TRUE)
  
  rsemGlmnet <- elasticNet(lavaanModel = modelFit, 
                           weights =  weights, 
                           alphas = 1, 
                           lambdas = lambdas,
                           method = "glmnet",
                           control = controlGlmnet(verbose = 0, 
                                                   startingValues = "est")
  )
  testthat::expect_equal(all(abs(rsemGlmnet@parameters[,rsemGlmnet@regularized] - 
                                   lslxParameter[,regularized]) < .002), TRUE)
  plot(rsemGlmnet)
  coef(rsemGlmnet)
  coef(rsemGlmnet, alpha = 1, lambda = lambdas[1])
  coef(rsemGlmnet, criterion = "AIC")
  coef(rsemGlmnet, criterion = "BIC")
  
  testthat::expect_equal(all(round(rsemGlmnet@fits$m2LL - (regsemM2LL))==0), TRUE)
  penalty <- N*lambdas*apply(regsemPars[,names(lavaanParameters)],1, function(x) sum(abs(x)*weights))
  testthat::expect_equal(all(round(rsemGlmnet@fits$regM2LL - (regsemM2LL+penalty))==0), TRUE)
  
  ## Test approximated cross-validation
  warning("Not testing approximate cross-validation.")
  cv <- aCV4regularizedSEM(regularizedSEM = rsem, k = N)
  # coef(cv)
  # coef(cv, alpha = 1, lambda = lambdas[1])
  # plot(cv)
  # 
  # # set automatic lambda:
  # rsem2 <- regularizeSEM(lavaanModel = modelFit, 
  #                        regularizedParameterLabels = regularizedLavaan,
  #                        penalty = "adaptiveLasso", 
  #                        lambdas = NULL,
  #                        nLambdas = 10)
  # testthat::expect_equal(all(apply(rsem2@parameters[,regularizedLavaan] == 0,2,sum) > 0), TRUE)
})
