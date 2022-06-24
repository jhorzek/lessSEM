test_that("testing general purpose ista lasso", {
  library(linr)
  set.seed(123)
  n <- 100
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    x4 = rnorm(n),
    x5 = rnorm(n)
  )
  df$y <- 0*df$x1 + .2*df$x2 + .3*df$x3 + .4*df$x4 + .5*df$x5 + rnorm(n,0,.3) 
  lmFit <- lm(y ~ 0+., df)
  
  # manual
  X <- as.matrix(df[,grepl("x", colnames(df))])
  listElements <- list(X = X, 
                       y = df$y)
  fitFunction <- function(b, parameterLabels, listElements){
    pred <- listElements$X%*%t(b)
    sse <- sum((listElements$y - pred)^2)
    return(sse)
  }
  listElements$fitFunction <- fitFunction
  
  gradientFunction <- function(b, parameterLabels, listElements){
    grad <- numDeriv::grad(func = listElements$fitFunction, 
                           x = b, 
                           listElements = listElements)
    names(grad) <- names(b)
    return(grad)
  }
  
  # initialize
  b <- rep(0, 5)
  names(b) <- paste0("b",1:5)
  fitFunction(matrix(b,nrow = 1), names(b), listElements)
  gradientFunction(matrix(b,nrow = 1), names(b), listElements)
  
  weights <- b
  weights[paste0("b",1:3)] <- 1
  
  control <- list(
    L0 = .1,
    eta = 2,
    maxIterOut = 10,
    maxIterIn = 20,
    breakOuter = .00000001,
    convCritInner = 1,
    sigma = .01,
    stepSizeInheritance = 3,
    verbose = 10
  )
  
  IL <- new(istaEnetGeneralPurpose, 
            weights, 
            control)
  
  lassoResult <- IL$optimize(
    fitFunction,
    gradientFunction,
    listElements,
    b,
    0,
    1
  )
  
  lassoResult$fit - sum((lmFit$residuals)^2)
  
  lassoResult$rawParameters
  
  testthat::expect_equal(all(round(lassoResult$rawParameters -
                                     coef(lmFit),
                                   4) == 0),
                         TRUE)
  
  control <- list(
    L0 = .1,
    eta = 2,
    maxIterOut = 1000,
    maxIterIn = 1000,
    breakOuter = .00000001,
    convCritInner = 1,
    sigma = .01,
    stepSizeInheritance = 0,
    verbose = 10
  )
  
  IL <- new(istaEnetGeneralPurpose, 
            weights, 
            control)
  
  lassoResult <- IL$optimize(
    fitFunction,
    gradientFunction,
    listElements,
    b,
    500,
    1
  )
  lassoResult$fit
  lassoResult$rawParameters
  coef(lmFit)
  testthat::expect_equal(all(round(lassoResult$rawParameters[c("b4", "b5")] -
                                     coef(lm(y ~ 0+x4+x5, df)),
                                   4) == 0),
                         TRUE)
})
