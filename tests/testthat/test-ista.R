test_that("testing gradients", {
  library(lavaan)
  library(aCV4SEM)
  
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
  
  SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = model, whichPars = "start")
  SEM <- aCV4SEM:::fit(SEM)
  SEM$getGradients(T)
  eigen(SEM$impliedCovariance)$values
  
  startingValues <- aCV4SEM:::getParameters(SEM, TRUE)
  labels <- SEM$getParameterLabels()
  
  weights <- startingValues
  weights[] <- 0
  weights[labels %in% c("a", "b")] <- 1
  
  lassoPar <- aCV4SEM:::istaLasso(SEM, 
                                  lambda = 0, 
                                  weights, 
                                  maxIterOut = 100000)
  SEM$fit() - -2*logLik(model)
  
  testthat::expect_equal(all(round(lassoPar - 
                                     aCV4SEM::getLavaanParameters(model)[names(lassoPar)],
                             1) == 0),
  TRUE)
})
