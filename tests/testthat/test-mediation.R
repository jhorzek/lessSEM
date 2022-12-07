test_that("testing mediation model", {
  testthat::skip_on_cran()
  library(lavaan) 
  library(lessSEM)
  N <- 100
  
  x_0 <- stats::rnorm(n = N, mean = 0, sd = 1)
  x_1 <- stats::rnorm(N, mean = 0, sd = 1)
  
  mediator <- .4*x_0 + .3*x_1 + stats::rnorm(n = N, mean = 0, sd = sqrt(.3))
  y <- .45*mediator + stats::rnorm(n = N, mean = 0, sd = sqrt(.3))
  
  datensatz <- data.frame(x_0 = x_0,
                          x_1 = x_1,
                          mediator = mediator,
                          y = y)
  
  lavaanModell <- "mediator ~ x_0 + x_1
  y ~ mediator 
x_0 ~~ x_0
x_1 ~~ x_1
mediator ~~ mediator
y ~~ y
"
  
  fitModel <- lavaan(model = lavaanModell, data = datensatz)
  
  regsem <- lasso(lavaanModel = fitModel, 
                  regularized = c("mediator~x_1"),
                  lambdas = 0)
  
  testthat::expect_equal(abs(regsem@fits$m2LL[1] - (-2*logLik(fitModel))) < 1e-3,TRUE)
  
})
