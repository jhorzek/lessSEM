test_that("testing multi-group SEM", {
  testthat::skip_on_cran()
  library(lessSEM)
  set.seed(123)
  library(lessSEM)
  
  # Identical to regsem, lessSEM builds on the lavaan
  # package for model specification. The first step
  # therefore is to implement the model in lavaan.
  
  dataset <- simulateExampleData(N = 200)
  
  lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"
  
  lavaanModel1 <- lavaan::sem(lavaanSyntax,
                              data = dataset[1:100,],
                              meanstructure = TRUE,
                              std.lv = TRUE)
  lavaanModel2 <- lavaan::sem(lavaanSyntax,
                              data = dataset[101:200,],
                              meanstructure = TRUE,
                              std.lv = TRUE)
  
  regsem <- bfgs(
    # pass the fitted lavaan model
    lavaanModel = c(lavaanModel1, lavaanModel2))
  regsem@parameters
  
  lavaanModelCombined <- lavaan::sem(lavaanSyntax,
                                     data = dataset,
                                     meanstructure = TRUE,
                                     std.lv = TRUE)
  regsem@fits$m2LL - -2*logLik(lavaanModelCombined)
  
  testthat::expect_equal(abs(regsem@fits$m2LL - -2*logLik(lavaanModelCombined)) > .0001,
                         FALSE)
})
