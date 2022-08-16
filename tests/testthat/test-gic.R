test_that("testing gic", {
  set.seed(1234)
  library(lessSEM)

  dataset <- simulateExampleData()
  
  lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"
  
  lavaanModel <- lavaan::sem(lavaanSyntax,
                             data = dataset,
                             meanstructure = TRUE,
                             std.lv = TRUE)
  
  regsem <- lasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    nLambdas = 50)
  
  AICs <- AIC(regsem)$AIC
  nPar <- AIC(regsem)$nonZeroParameters
  
  regsemSmooth <- smoothLasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    lambdas = regsem@fits$lambda,
    epsilon = 1e-8,
    tau = 0)
  
  GAIC <- GIC(regularizedSEM = regsemSmooth, k = 2)
  
  testthat::expect_equal(all(abs(nPar - GAIC$df) < 2), TRUE)
})
