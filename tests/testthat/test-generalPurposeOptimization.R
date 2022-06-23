test_that("testing general purpose ista lasso", {
  library(linr)
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  x5 <- rnorm(n)
  
  y <- 0*x1+.2*x2+.3*x3 + .4*x4 + .5*x5 + rnorm(n,0,.3) 
  df <- data.frame(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x5 = x5,
    y = y
  )
  lmFit <- lm(y ~ 0+., df)
  
  # manual
  X <- as.matrix(df[,grepl("x", colnames(df))])
  listElements <- list(X = X, 
                      y = y)
  fitFunction <- function(b, listElements){
    pred <- listElements$X%*%b
    sse <- sum((listElements$y - pred)^2)
    return(sse)
  }
  listElements$fitFunction <- fitFunction
  
  gradientFunction <- function(b, listElements){
    grad <- numDeriv::grad(func = listElements$fitFunction, 
                           x = b, 
                           listElements = listElements)
    names(grad) <- names(b)
    return(grad)
  }
  
  # initialize
  b <- rep(0, 5)
  names(b) <- paste0("b",1:5)
  fitFunction(b, listElements)
  gradientFunction(b, listElements)
  
  weights <- b
  weights[] <- 1
  
  control <- list(
    i_k = 2,
    L0 = 1,
    eta = 2,
    maxIterOut = 50,
    maxIterIn = 500,
    breakOuter = .00000001
  )
  
  IL <- new(istaLASSOGeneralPurpose, weights, control)
  
  lassoResult <- IL$optimize(
    fitFunction,
    gradientFunction,
    listElements,
    b,
    0
  )
  
  lassoResult$fit - sum((lmFit$residuals)^2)
  
  lassoResult$rawParameters
  
  testthat::expect_equal(all(round(lassoResult$rawParameters -
                                     coef(lmFit),
                                   4) == 0),
                         TRUE)
  
  lassoResult <- IL$optimize(
    fitFunction,
    gradientFunction,
    listElements,
    b,
    n*.6
  )
  lassoResult$fit
  lassoResult$rawParameters
  coef(lmFit)
})
