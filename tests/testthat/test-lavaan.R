test_that("testing lavaan", {
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
  
  model <- sem(model1, data = PoliticalDemocracy, meanstructure = TRUE)
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model)
  show(SEM)
  logLik(SEM)
  testthat::expect_equal(round(AIC(SEM) - AIC(model),5),0)
  testthat::expect_equal(round(BIC(SEM) - BIC(model),5),0)
  
  testthat::expect_equal(round(SEM$m2LL - (-2*as.numeric(logLik(model))),4),0)
  testthat::expect_equal(round(model@Fit@test$standard$stat - lessSEM:::.likelihoodRatioFit(par = lessSEM:::.getParameters(SEM), SEM),4)[1,1],0)

  # test missing data
  dat <- PoliticalDemocracy
  
  dat[1:3,4] <- NA
  
  dat[12,5] <- NA
  
  model <- sem(model1, data = dat, meanstructure = TRUE, missing = "ML")
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model)
  lessSEM:::.fit(SEM)
  
  testthat::expect_equal(round(SEM$m2LL - (-2*as.numeric(logLik(model))),4),0)

})
