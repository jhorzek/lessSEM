test_that("testing standard errors", {
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  model1 <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8

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
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model)
  SEM <- lessSEM:::.fit(SEM)
  
  se <- lessSEM:::.standardErrors(SEM = SEM, raw = FALSE)
  
  seLavaan <- parameterEstimates(model)[,"se"]
  seLavaan <- seLavaan[seLavaan != 0]
  
  testthat::expect_equal(all(abs(se - seLavaan) < .01), TRUE)
})
