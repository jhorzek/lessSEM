test_that("testing penalty interfaces", {
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  
  # Identical to regsem, lessSEM builds on the lavaan
  # package for model specification. The first step
  # therefore is to implement the model in lavaan.
  set.seed(123)
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
  
  regsem <- try(lasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    method = "ista",
    control = controlIsta())
  )
  
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(lasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    method = "ista",
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)

  
  regsem <- try(lasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    method = "glmnet",
    control = controlGlmnet())
  )
  
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(lasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    nLambdas = 5,
    method = "glmnet",
    control = controlGlmnet())
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(ridge(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    method = "ista",
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(ridge(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    method = "glmnet",
    control = controlGlmnet())
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(adaptiveLasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    method = "ista",
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  
  regsem <-  try(adaptiveLasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    nLambdas = 5,
    method = "ista",
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  
  regsem <-  try(adaptiveLasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    method = "glmnet",
    control = controlGlmnet())
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(adaptiveLasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    nLambdas = 5,
    method = "glmnet",
    control = controlGlmnet())
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  
  regsem <-  try(elasticNet(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    alphas = seq(0,1,length.out = 3),
    method = "ista",
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(elasticNet(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    alphas = seq(0,1,length.out = 3),
    method = "glmnet",
    control = controlGlmnet())
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  
  regsem <-  try(cappedL1(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    thetas = seq(0.1,1,length.out = 3), 
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(lsp(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    thetas = seq(0.1,1,length.out = 3), 
    control = controlIsta())
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(mcp(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    thetas = seq(0.1,1,length.out = 3), 
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
  
  regsem <-  try(scad(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    lambdas = seq(0,1,length.out = 5),
    thetas = seq(2.1,3,length.out = 3), 
    control = controlIsta()) 
  )
  testthat::expect_equal(is(regsem, "try-error"), FALSE)
}
)
