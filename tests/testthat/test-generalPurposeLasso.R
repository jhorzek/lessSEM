test_that("testing general purpose ista lasso", {
  library(lessSEM)
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
    accelerate = TRUE,
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
  
  unregularized <- IL$optimize(
    b,
    fitFunction,
    gradientFunction,
    listElements,
    0,
    1
  )
  
  unregularized$fit - sum((lmFit$residuals)^2)
  
  unregularized$rawParameters
  
  testthat::expect_equal(all(round(unregularized$rawParameters -
                                     coef(lmFit),
                                   4) == 0),
                         TRUE)
  
  control <- list(
    L0 = .1,
    eta = 2,
    accelerate = TRUE,
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
    b,
    fitFunction,
    gradientFunction,
    listElements,
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
  
  control <- list(
    initialHessian = diag(1,length(b)),
    stepSize = .9,
    sigma = 0,
    gamma = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 1000,
    breakOuter = 1e-6,
    breakInner = 1e-5,
    convergenceCriterion = 0, 
    verbose = 0
  )
  IL <- new(glmnetEnetGeneralPurpose, 
            weights, 
            control)
  
  lassoResult <- IL$optimize(
    b,
    fitFunction,
    gradientFunction,
    listElements,
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
