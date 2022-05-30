test_that("testing optimization with PoliticalDemocracy", {
  library(lavaan)
  library(aCV4SEM)
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
  
  regsem <- regularizeSEM(lavaanModel = model, 
                          regularizedParameterLabels = c("a", "b", "c"),
                          penalty = "lasso",
                          lambdas = 0, 
                          control = controlGLMNET(startingValues = "start"))
  testthat::expect_equal(abs(regsem@fits$m2LL[1] - (-2*logLik(model))) < 1e-3,TRUE)
  
})
