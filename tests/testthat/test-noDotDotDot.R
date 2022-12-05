test_that("testing noDotDotDot", {
  testthat::skip_on_cran()
  library(lessSEM)
  
  # test function that only takes one argument
  fn <- function(params){
    return(params)
  }
  
  out <- try(
    lessSEM:::.noDotDotDot(fn = fn, fnName = "fn")
  )
  
  testthat::expect_equal(is(out, "list"), TRUE)
  
  # try calling the function
  
  testthat::expect_equal(out$fnUser(1, "a", additionalArguments = out$additionalArguments) == 1, 
                         c(a = TRUE))
  
  # test function that takes additional arguments
  fn <- function(params, x, y, z){
    return(params)
  }
  
  out <- try(
    lessSEM:::.noDotDotDot(fn = fn, fnName = "fn", x = 1, y = "a", z = 1)
  )
  
  testthat::expect_equal(is(out, "list"), TRUE)
  
  # try calling the function
  
  testthat::expect_equal(out$fnUser(1, "a", additionalArguments = out$additionalArguments) == 1, 
                         c(a = TRUE))
  
})
