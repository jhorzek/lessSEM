test_that("testing regsem", {
  library(regsem)
  library(aCV4SEM)
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = FALSE)
  # increase to > 25
  cv.out = cv_regsem(outt,type="lasso", pars_pen=c(1:2,6:8),
                     n.lambda=10,jump=0.01, round = 10)
  
  parameters <- cvregsem2LavaanParameters(cvregsemModel = cv.out, lavaanModel = outt)
  
  CFA <- SEMFromLavaan(lavaanModel = outt)

  chisquares <- rep(NA, nrow(parameters))
  for(p in 1:nrow(parameters)){
    chisquares[p] <- likelihoodRatioFit(par = parameters[p,], SEM = CFA, raw = FALSE)
  }
  
  testthat::expect_equal(sum(round(chisquares - cv.out$fits[,"chisq"],4)), 0)
  
  
})
