test_that("meanstructure works", {
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
  modelFitNoMeans = cfa(modelSyntax, y, meanstructure = FALSE)
  modelFitWithMeans <- cfa(modelSyntax, y, meanstructure = TRUE)
  
  lambdas <- seq(0,1,length.out = 5)
  
  rsemGlmnetNoMeans <- lasso(lavaanModel = modelFitNoMeans, 
                             regularized = paste0("f=~y",6:ncol(y)), 
                             lambdas = lambdas,
                             method = "glmnet",
                             control = controlGlmnet()
  )
  
  rsemGlmnetWithMeans <- lasso(lavaanModel = modelFitNoMeans, 
                               regularized = paste0("f=~y",6:ncol(y)), 
                               lambdas = lambdas,
                               method = "glmnet",
                               control = controlGlmnet(),
                               modifyModel = modifyModel(addMeans = TRUE)
  )
  
  
  testthat::expect_equal(all(round(rsemGlmnetWithMeans@parameters[,colnames(rsemGlmnetNoMeans@parameters)] - 
                                     rsemGlmnetNoMeans@parameters,2)==0), TRUE)
  
  rsemIstaNoMeans <- lasso(lavaanModel = modelFitNoMeans, 
                           regularized = paste0("f=~y",6:ncol(y)), 
                           lambdas = lambdas,
                           method = "ista",
                           control = controlIsta()
  )
  
  rsemIstaWithMeans <- lasso(lavaanModel = modelFitNoMeans, 
                             regularized = paste0("f=~y",6:ncol(y)), 
                             lambdas = lambdas,
                             method = "ista",
                             control = controlIsta(),
                             modifyModel = modifyModel(addMeans = TRUE)
  )
  
  testthat::expect_equal(all(round(rsemIstaWithMeans@parameters[,colnames(rsemIstaNoMeans@parameters)] - 
                                     rsemIstaNoMeans@parameters,4)==0), TRUE)
  
  modelSyntax2 <- paste0(modelSyntax, "\n", "f ~ 1", "\n", paste0(paste0(colnames(y), " ~ 0"), collapse = "; "))
  
  modelFit2WithMeans <- cfa(modelSyntax2, y)
  
  rsemGlmnet2WithMeans <- lasso(lavaanModel = modelFit2WithMeans, 
                                regularized = paste0("f=~y",6:ncol(y)), 
                                lambdas = 0,
                                method = "glmnet",
                                control = controlGlmnet()
  )
  
  testthat::expect_equal(all(round(rsemGlmnet2WithMeans@parameters[,names(getLavaanParameters(modelFit2WithMeans))] -
                                     getLavaanParameters(modelFit2WithMeans), 3) == 0),
                         TRUE)
  
  rsemIsta2WithMeans <- lasso(lavaanModel = modelFit2WithMeans, 
                              regularized = paste0("f=~y",6:ncol(y)), 
                              lambdas = 0,
                              method = "ista",
                              control = controlIsta()
  )
  
  testthat::expect_equal(all(round(rsemIsta2WithMeans@parameters[,names(getLavaanParameters(modelFit2WithMeans))] -
                                     getLavaanParameters(modelFit2WithMeans), 3) == 0),
                         TRUE)
  
})
