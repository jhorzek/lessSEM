test_that("optimization works", {
  
  library(regsem)
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
  
  CFA <- SEMFromLavaan(lavaanModel = modelFit, rawData = y)
  CFA <- fit(CFA)
  correctParameters <- getParameters(CFA, raw = FALSE)
  
  pars <- getParameters(CFA, raw  = TRUE)
  pars <- pars + rnorm(length(pars), 0, .3)
  CFA <- setParameters(CFA, names(pars), pars, raw = TRUE)
  opt <- optimizeSEM(SEM = CFA, raw = TRUE)
  testthat::expect_equal(any(abs(getParameters(opt$SEM, raw = FALSE) - correctParameters) > .001), FALSE)
  
  # test lasso penalty
  
  lambda <- .5
  pars_pen = 5:24
  
  regsem_fit <- regsem(model = modelFit, 
                       lambda = lambda, 
                       pars_pen = pars_pen, 
                       gradFun = "ram")
  
  pars <- regsem2LavaanParameters(regsemModel = regsem_fit, lavaanModel = modelFit)
  regularizedParameters <- names(pars)[pars_pen]
  
  optLasso <- optimizeRegularizedSEM(lavaanModel = modelFit,
                                     regularizedParameterLabels = regularizedParameters, 
                                     lambda = lambda, 
                                     penalty = "lasso")
  testthat::expect_equal(any(abs(optLasso$parameters[,names(pars)] - pars) > .1),
                         FALSE)
  
  regsem_cvFit <- cv_regsem(model = modelFit, 
                            pars_pen = pars_pen, 
                            gradFun = "ram")
  pars_cvFit <- cvregsem2LavaanParameters(cvregsemModel = regsem_cvFit, lavaanModel = modelFit)
  
  optLasso2 <- optimizeRegularizedSEM(lavaanModel = modelFit,
                                     regularizedParameterLabels = regularizedParameters, 
                                     lambda = regsem_cvFit$fits[,"lambda"], 
                                     penalty = "lasso")
  testthat::expect_equal(any(abs(optLasso2$parameters[,names(pars)] - pars_cvFit[,names(pars)]) > .1),
                         FALSE)

})
