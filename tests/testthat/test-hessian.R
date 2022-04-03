test_that("hessian works", {
  library(lavaan)
  
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
  N <- nrow(PoliticalDemocracy)
  PoliticalDemocracy[1:4,2] <- NA
  model <- sem(model1, data = PoliticalDemocracy, meanstructure = TRUE, missing = "ML")
  parEst <- parameterEstimates(model, remove.nonfree = TRUE)
  parLabels <- paste0(parEst$lhs, parEst$op, parEst$rhs)
  parLabels[parEst$label != ""] <- parEst$label[parEst$label != ""]
  
  SEM <- SEMFromLavaan(model = model, rawData = PoliticalDemocracy)
  SEM <- fit(SEM)
  
  grad <- computeAnalyticGradients(par = getParameters(SEM), SEM = SEM)
  
  hessian <- computeHessianFromAnalytic(par = getParameters(SEM), SEM = SEM, eps = 1e-8)
  
  hessian2 <- numDeriv::hessian(func = minus2LogLikelihood, x = getParameters(SEM), SEM = SEM, raw = FALSE)
  
  testthat::expect_equal(abs(sum(hessian - hessian2)) < .001, TRUE)
  
})
