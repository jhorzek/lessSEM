test_that("testing automatic standardization for cross-validation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("regsem")
  library(regsem)
  library(lessSEM)
  set.seed(123)
  N <- 50
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
  
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  rsem <- lessSEM::lasso(lavaanModel = modelFit, 
                         regularized = regularizedLavaan,
                         nLambdas = 5)
  
  lambdas <- rsem@fits$lambda
  
  # test standardization
  cv <- cvLasso(lavaanModel = modelFit, 
                regularized = regularizedLavaan,
                lambdas = lambdas, 
                returnSubsetParameters = TRUE,
                standardize = TRUE,
                k = 2)
  
  subsets <- cv@subsets
  pars <- cv@subsetParameters
  cvfits <- cv@cvfits
  
  parameterLabels <- cv@parameterLabels  
  
  # lavaan tends to sort the data differently
  ySorted <- try(lavaan::lavInspect(modelFit, 
                                    "data"))

  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = modelFit)
  
  for(ro in 1:nrow(pars)){
    
    trainSet <- ySorted[!subsets[,pars$trainSet[ro]],]
    testSet <- ySorted[subsets[,pars$trainSet[ro]],]
    
    means <- apply(trainSet,2,mean)
    standardDeviations <- apply(trainSet,2,sd)
    
    testSet <- lessSEM::cvScaler(testSet = testSet, 
                                 means = means, 
                                 standardDeviations = standardDeviations)
    
    SEM <- lessSEM:::.setParameters(SEM = SEM, 
                                    labels = parameterLabels, 
                                    values = unlist(pars[ro, parameterLabels]),
                                    raw = FALSE)
    SEM$fit()
    
    m2LL <- -2*sum(mvtnorm::dmvnorm(
      x = testSet, 
      mean = SEM$impliedMeans,
      sigma = SEM$impliedCovariance,
      log = TRUE
    ))
    sel <- cvfits$lambda == pars$lambda[ro]
    if(sum(sel) != 1) stop("Error when selecting cv target")
    cvfits[sel,"cvfit"] <- cvfits[sel,"cvfit"] - m2LL
    
  }
  
  testthat::expect_equal(all(abs(cvfits$cvfit)< 1e-6), TRUE)
  
})
