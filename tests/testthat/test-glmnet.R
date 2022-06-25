test_that("testing ista-lasso", {
  library(lavaan)
  library(linr)
  
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
  
  SEM <- linr:::SEMFromLavaan(lavaanModel = model, whichPars = "start")
  SEM <- linr:::fit(SEM)
  SEM$getGradients(T)
  eigen(SEM$impliedCovariance)$values
  
  startingValues <- linr:::getParameters(SEM, TRUE)
  labels <- SEM$getParameterLabels()
  
  weights <- startingValues
  weights[] <- 0
  weights[labels %in% c("a", "b")] <- 1
  
  control <- list(
    initialHessian = diag(1,length(startingValues)),
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
  
  GNet <- new(glmnetEnet, weights, control)
  unregularized <- GNet$optimize(
    SEM,
    startingValues,
    0,
    1
  )
  
  unregularized$fit - -2*logLik(model)
  
  unregularized$rawParameters
  
  testthat::expect_equal(all(round(linr:::getParameters(SEM, raw = FALSE) - 
                                     linr::getLavaanParameters(model)[names(linr:::getParameters(SEM, raw = FALSE))],
                                   1) == 0),
                         TRUE)
  
  # lasso
  lassoResult <- GNet$optimize(
    SEM,
    startingValues,
    1.5,
    1
  )
  lassoResult$rawParameters
  testthat::expect_equal(any(lassoResult$rawParameters == 0),
                         TRUE)
  
  # enet
  enetResult <- GNet$optimize(
    SEM,
    startingValues,
    2.7,
    .4
  )
  enetResult$rawParameters
  
  control <- list(
    L0 = .1,
    eta = 2,
    maxIterOut = 10000,
    maxIterIn = 1000,
    breakOuter = .000000001,
    convCritInner = 1,
    sigma = .1,
    stepSizeInheritance = 3,
    verbose = 0
  )
  
  IN <- new(istaEnet, weights, control)
  lassoResultIsta <- IN$optimize(
    SEM,
    startingValues,
    1.5,
    1
  )
  lassoResultIsta$rawParameters
  round(lassoResult$rawParameters -
  lassoResultIsta$rawParameters,3)
  round(lassoResult$fit -
          lassoResultIsta$fit,7)
  
  enetResultIsta <- IN$optimize(
    SEM,
    startingValues,
    2.7,
    .4
  )
  enetResultIsta$rawParameters
  
  enetResult$rawParameters - enetResultIsta$rawParameters
  enetResult$fit - enetResultIsta$fit
})
