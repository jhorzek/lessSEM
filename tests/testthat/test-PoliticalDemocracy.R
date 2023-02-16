test_that("testing optimization with PoliticalDemocracy", {
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  modelSyntax <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
  
  model <- sem(modelSyntax, 
               data = PoliticalDemocracy, 
               meanstructure = TRUE)
  
  modelNoFit <- sem(modelSyntax, 
                    data = PoliticalDemocracy, 
                    meanstructure = TRUE,
                    do.fit = FALSE)
  
  regsem <- lasso(lavaanModel = modelNoFit, 
                  regularized = c("a", "b", "c"),
                  lambdas = 0,
                  method = "ista",
                  control = controlIsta(startingValues = "start"))
  testthat::expect_equal(abs(regsem@fits$m2LL[1] - (-2*logLik(model))) < 1e-3,TRUE)
  
  regsem <- lasso(lavaanModel = modelNoFit, 
                  regularized = c("a", "b", "c"),
                  lambdas = 0,
                  control = controlGlmnet(startingValues = "start", initialHessian = "compute", verbose = 1, maxIterOut = 100))
  testthat::expect_equal(abs(regsem@fits$m2LL[1] - (-2*logLik(model))) < 1e-3,TRUE)
  
  regsem <- bfgs(lavaanModel = modelNoFit,
                 control = controlBFGS(startingValues = "start"))
  testthat::expect_equal(abs(regsem@fits$m2LL[1] - (-2*logLik(model))) < 1e-3,TRUE)
  
})
