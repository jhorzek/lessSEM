test_that("lasso works", {
  library(lslx)
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
  fitLslx$fit(penalty_method = "lasso",lambda_grid = lambdas, loss = "ml")
  
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
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  rsem <- regularizeSEM(lavaanModel = modelFit, 
                        regularizedParameterLabels = regularizedLavaan,
                        penalty = "lasso", 
                        lambdas = lambdas)
  testthat::expect_equal(all(round(wideResults(rsem)[,regularizedLavaan] - lslxParameter[,regularized],3)==0), TRUE)
  plot(rsem)
  coef(rsem)
  coef(rsem, alpha = 1, lambda = .1)
  
  ## Test approximated cross-validation
  
  cv <- aCV4regularizedSEM(regularizedSEM = rsem, k = N)
  coef(cv)
  coef(cv, alpha = 1, lambda = .1)
  plot(cv)
})
