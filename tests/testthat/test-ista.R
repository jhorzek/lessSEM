test_that("testing ista-lasso", {
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  
  model1 <- ' 
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
  
  model <- sem(model1, 
               data = PoliticalDemocracy, 
               meanstructure = TRUE)
  scores <- -2*lavScores(model)
  gradients <- apply(scores,2,sum)
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model, whichPars = "start")
  SEM <- lessSEM:::.fit(SEM)
  SEM$getGradients(T)
  eigen(SEM$impliedCovariance)$values
  
  startingValues <- lessSEM:::.getParameters(SEM, TRUE)
  labels <- SEM$getParameterLabels()
  
  weights <- startingValues
  weights[] <- 0
  weights[labels %in% c("a", "b")] <- 1
  
  control <- list(
    L0 = .1,
    eta = 2,
    accelerate = TRUE,
    maxIterOut = 2000,
    maxIterIn = 1000,
    breakOuter = .00000001,
    convCritInner = 1,
    sigma = .01,
    stepSizeInheritance = 3,
    verbose = 0
  )
  
  IL <- new(istaEnetSEM, weights, control)
  lassoResult <- IL$optimize(
    startingValues,
    SEM,
    0,
    1
  )
  
  lassoResult$fit - -2*logLik(model)
  
  lassoResult$rawParameters
  
  testthat::expect_equal(all(round(lessSEM:::.getParameters(SEM, raw = FALSE) - 
                                     lessSEM::getLavaanParameters(model)[names(lessSEM:::.getParameters(SEM, raw = FALSE))],
                                   1) == 0),
                         TRUE)
  
  # lasso
  lassoResult <- IL$optimize(
    startingValues,
    SEM,
    100,
    1
  )
  lassoResult$rawParameters
  testthat::expect_equal(any(lassoResult$rawParameters == 0),
                         TRUE)
  
  control <- list(
    L0 = .1,
    eta = 2,
    accelerate = FALSE,
    maxIterOut = 10000,
    maxIterIn = 1000,
    breakOuter = .00000001,
    convCritInner = 1,
    sigma = .1,
    stepSizeInheritance = 3,
    verbose = 0
  )
  
  IL <- new(istaEnetSEM, weights, control)
  
  # enet
  lassoResult <- IL$optimize(
    startingValues,
    SEM,
    100,
    .4
  )
  lassoResult$rawParameters
})
