test_that("testing ista-lasso", {
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  
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
  
  H <- diag(1000, length(startingValues))
  control <- list(
    initialHessian = H,
    stepSize = .9,
    sigma = 0,
    gamma = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 1000,
    breakOuter = 1e-10,
    breakInner = 1e-10,
    convergenceCriterion = 0, 
    verbose = 0
  )
  
  GNet <- new(glmnetEnet, weights, control)
  unregularized <- GNet$optimize(
    startingValues,
    SEM,
    0,
    1
  )
  control$H <- unregularized$Hessian
  unregularized$fit - -2*logLik(model)
  
  unregularized$rawParameters
  unregularizedParam <- lessSEM:::.getParameters(SEM, raw = FALSE)
  
  testthat::expect_equal(all(round( unregularizedParam - 
                                     getLavaanParameters(model)[names(lessSEM:::.getParameters(SEM, raw = FALSE))],
                                   1) == 0),
                         TRUE)
  
  # lasso
  lambda_lasso <- 5
  lassoResult <- GNet$optimize(
    startingValues,
    SEM,
    lambda_lasso,
    1
  )
  lassoResult$rawParameters
  lassoParam <- lessSEM:::.getParameters(SEM, raw = FALSE)
  testthat::expect_equal(any(lassoResult$rawParameters == 0),
                         TRUE)
  
  # enet
  lambda_enet <- 1.5
  enetResult <- GNet$optimize(
    startingValues,
    SEM,
    lambda_enet,
    .4
  )
  enetResult$rawParameters
  enetParam <- lessSEM:::.getParameters(SEM, raw = FALSE)
  
  control <- list(
    L0 = .1,
    eta = 2,
    accelerate = TRUE,
    maxIterOut = 10000,
    maxIterIn = 1000,
    breakOuter = .00000000001,
    convCritInner = 1,
    sigma = .1,
    stepSizeInheritance = 3,
    verbose = 0
  )
  
  IN <- new(istaEnet, weights, control)
  lassoResultIsta <- IN$optimize(
    startingValues,
    SEM,
    lambda_lasso,
    1
  )
  lassoResultIsta$rawParameters
  lassoParamIsta <- lessSEM:::.getParameters(SEM, raw = FALSE)
  
  testthat::expect_equal(all(
    round(lassoParam -
            lassoParamIsta,2) == 0), TRUE)
  testthat::expect_equal(round(lassoResult$fit -
          lassoResultIsta$fit,3) == 0, TRUE)
  
  enetResultIsta <- IN$optimize(
    startingValues,
    SEM,
    lambda_enet,
    .4
  )
  enetResultIsta$rawParameters
  enetParamIsta <- lessSEM:::.getParameters(SEM, raw = FALSE)
  
  testthat::expect_equal(all(
    round(
      enetParam - enetParamIsta,2) == 0), TRUE)
  testthat::expect_equal(abs(enetResult$fit - enetResultIsta$fit) < .01, TRUE)
  
  N <- nrow(PoliticalDemocracy)
  
  SEM <- lessSEM:::.setParameters(SEM, names(enetResultIsta$rawParameters), 
                             enetResultIsta$rawParameters, raw = TRUE)
  testthat::expect_equal(abs(
    enetResultIsta$fit -
                           (SEM$fit() + 
    N*lambda_enet*.4* sum(abs(enetResultIsta$rawParameters*weights)) +
    N*lambda_enet*(1-.4) * sum((enetResultIsta$rawParameters*weights)^2))) < 1e-8, TRUE)
  
  
  SEM <- lessSEM:::.setParameters(SEM, names(enetResult$rawParameters), enetResult$rawParameters, raw = TRUE)
  
  testthat::expect_equal(
    enetResult$fit - (
  SEM$fit() + 
    N*lambda_enet*.4*sum(abs(enetResult$rawParameters*weights)) +
    N*lambda_enet*(1-.4) * sum((enetResult$rawParameters*weights)^2)
    ) == 0, TRUE)

})

