test_that("regsem works", {
  library(regsem)
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = FALSE)
  # increase to > 25
  cv.out = cv_regsem(outt,type="lasso", pars_pen=c(1:2,6:8),
                     n.lambda=10,jump=0.01, round = 10)
  
  CFA <- SEMFromLavaan(outt, rawData = HS)
  CFA <- fit(CFA)
  
  chisquares <- rep(NA, nrow(cv.out$parameters))
  for(p in 1:nrow(cv.out$parameters)){
    CFA2 <- setParameters(CFA, labels = names(cv.out$parameters[p,]), 
                         values = cv.out$parameters[p,],
                         labelsFrom = "regsem", 
                         raw = FALSE)
    chisquares[p] <- likelihoodRatioFit(par = getParameters(CFA2), SEM = CFA, raw = FALSE)
  }
  
  testthat::expect_equal(sum(round(chisquares - cv.out$fits[,"chisq"],4)), 0)
  
  model_out = regsem(outt,
                     type="lasso", 
                     pars_pen=c(1:2,6:8),
                     lambda = .05,
                     gradFun = "ram",
                     round = 10)
  
  CFA2 <- setParameters(CFA, labels = names(model_out$out$pars), 
                        values = unlist(model_out$out$pars),
                        labelsFrom = "regsem", 
                        raw = FALSE)
  
  o2 <- optimizeLassoSEM(SEM = CFA, 
                         regularizedParameters = names(getParameters(CFA))[model_out$pars_pen], 
                         lambda = .05)
  
  testthat::expect_equal(sum(abs(getParameters(CFA2) - getParameters(o2$SEM))) < .01, TRUE)
  
  ## Test individual chisquares
  N <- nrow(HS)
  indchisL <- rep(NA, N)
  for(i in 1:N){
    indchisL[i] <- individualLikelihoodRatioFitRegsemLASSO(par = getParameters(o2$SEM), 
                                                           SEM = CFA, 
                                                           data = HS[i,], 
                                                           regularizedParameterLabels = names(model_out$coefficients[, model_out$pars_pen]), 
                                                           lambda_ = .05, 
                                                           raw = FALSE,
                                                           eps = .0001)
    
  }
  testthat::expect_equal(sum(abs(model_out$fit - sum(indchisL))) < .001, TRUE)
  
  
})
