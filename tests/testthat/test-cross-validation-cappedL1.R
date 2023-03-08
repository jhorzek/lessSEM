test_that("testing cross-validation for cappedL1", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("regsem")
  library(regsem)
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
  
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  
  lambdas <- seq(0,1,length.out = 3)
  thetas <- c(.3,1)

  ## Test cross-validation
  
  cv <- cvCappedL1(lavaanModel = modelFit, 
                   regularized = regularizedLavaan,
                   lambdas = lambdas,
                   thetas = thetas,
                   returnSubsetParameters = TRUE,
                   k = 2)
  
  testthat::expect_equal(ncol(cv@subsets), 2)
  
  coef(cv)
  plot(cv)
  
  subsets <- cv@subsets
  pars <- cv@subsetParameters
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = modelFit)
  
  parameterLabels <- cv@parameterLabels
  
  # lavaan tends to sort the data differently
  ySorted <- try(lavaan::lavInspect(modelFit, 
                                    "data"))
  
  for(ro in 1:nrow(pars)){
    
    trainSet <- pars$trainSet[ro]
    testSet <- subsets[,trainSet]
    
    SEM <- lessSEM:::.setParameters(SEM = SEM, 
                                  labels = parameterLabels, 
                                  values = unlist(pars[ro, parameterLabels]),
                                  raw = FALSE)
    SEM$fit()
    
    m2LL <- -2*sum(mvtnorm::dmvnorm(
      x = ySorted[testSet,], 
      mean = SEM$impliedMeans,
      sigma = SEM$impliedCovariance,
      log = TRUE
    ))
    sel <- cv@cvfits$lambda == pars$lambda[ro] & cv@cvfits$theta == pars$theta[ro]
    if(sum(sel) != 1) stop("Error when selecting cv target")
    cv@cvfits[sel,"cvfit"] <- cv@cvfits[sel,"cvfit"] - m2LL
    
  }
  
  testthat::expect_equal(all(abs(cv@cvfits$cvfit)< 1e-6), TRUE)
  
  # test subset parameters
  subset <- sample(1:2, 1)
  subsetPars <- pars[pars$trainSet == subset,]
  
  subsetCappedL1 <- cappedL1(lavaanModel = modelFit, 
                             regularized = regularizedLavaan,
                             lambdas = lambdas,
                             thetas = thetas,
                             modifyModel = modifyModel(dataSet = y[!subsets[,subset],]))
  
  testthat::expect_equal(all(abs(subsetCappedL1@parameters - subsetPars[,colnames(subsetCappedL1@parameters)])< 1e-3), TRUE)
  
}
)