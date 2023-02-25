test_that("unidentified works", {
  testthat::skip_on_cran()
  set.seed(123)
  library(lessSEM)
  
  sim <- "
  f1 =~ 1*y1 + .8*y2 + .7*y3
  f2 =~ 1*y4 + .8*y5 + .7*y6
  "
  data <- lavaan::simulateData(sim)
  
  model <- "
  f1 =~ 1*y1 + y2 + y3 + a*y4 + b*y5 + c*y6
  f2 =~ 1*y4 + d*y1 + e*y2 + f*y3 + y5 + y6
  "
  fit <- suppressWarnings(sem(model = model, 
                              data = data))
  
  reg <- lessSEM::lasso(fit,
                        regularized = c("a", "b", "c", "d", "e", "f"),
                        nLambdas = 10)
})
