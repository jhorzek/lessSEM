test_that("testing transformations", {
  library(lavaan)
  library(lessSEM)
  model1 <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + d*y6 + e*y7 + f*y8

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
  
  transformations <- "
  parameters: a,b,c,d,e,f,deltaA, deltaB, deltaC
  d = a + deltaA
  e = b + deltaB
  f = c + deltaC
  "
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model, 
                                  whichPars = "start",
                                  fit = FALSE, 
                                  addMeans = TRUE,
                                  transformations = transformations)
  show(SEM)
  logLik(SEM)
  
  params <- lessSEM:::.getParameters(SEM, raw = TRUE)
  testthat::expect_equal((params["d"] - params["a"] - params["deltaA"]) == 0,c("d" = TRUE))
  testthat::expect_equal((params["e"] - params["b"] - params["deltaB"]) == 0,c("e" = TRUE))
  testthat::expect_equal((params["f"] - params["c"] - params["deltaC"]) == 0,c("f" = TRUE))
  
  individualFit <- rep(NA, nrow(PoliticalDemocracy))
  for(i in 1:nrow(PoliticalDemocracy)){
    individualFit[i] <- lessSEM:::.individualMinus2LogLikelihood(par = lessSEM:::.getParameters(SEM), 
                                                                 SEM = SEM, 
                                                                 data = SEM$rawData[i,], 
                                                                 raw = FALSE)
  }
  
  testthat::expect_equal(round(sum(individualFit) - (-2*as.numeric(logLik(model))),4),0)
  
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
