test_that("testing approximate cross validation", {
  
  library(lavaan)
  library(regsem)
  library(aCV4SEM)
  #### Test lavaan ####
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = TRUE)
  
  aLOOCV <- aCV4lavaan(lavaanModel = outt, k = nrow(HS), raw = FALSE)
  
  exactLOOCV <- rep(NA, nrow(HS))
  for(i in 1:nrow(HS)){
    outt = cfa(mod, HS[-i,], meanstructure = TRUE)
    exactLOOCV[i] <- aCV4SEM:::computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                                     rawData = as.numeric(HS[i,]),
                                                     impliedMeans = outt@implied$mean[[1]], 
                                                     impliedCovariance = outt@implied$cov[[1]])
  }
  
  testthat::expect_equal(max(abs(aLOOCV$leaveOutFits -exactLOOCV)) < .1, TRUE)
  
  plot(exactLOOCV, exactLOOCV, type = "l",
       xlab = "exact loocv", ylab = "approximated loocv")
  points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
  
  #### Test regsem ####
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = TRUE)
  
  CFA <- aCV4SEM:::SEMFromLavaan(lavaanModel = outt)
  CFA <- aCV4SEM:::fit(CFA)
  
  lambda_ <- 0.05
  model_out = regsem(outt,
                     type="lasso", 
                     pars_pen=c(1:2,6:8),
                     lambda = lambda_,
                     gradFun = "ram",
                     round = 10)
  pars <- unlist(model_out$coefficients)
  
  # change labels to those of lavaan
  labels <- lavaan2regsemLabels(lavaanModel = outt)
  pars <- pars[labels$regsemLabels]
  names(pars) <- labels$lavaanLabels
  
  CFA <- aCV4SEM:::setParameters(CFA, 
                                 labels = names(pars), 
                                 values = pars, 
                                 raw = FALSE)
  fitDiff <- .5*(1/nrow(HS))*aCV4SEM:::likelihoodRatioFit(par = aCV4SEM:::getParameters(CFA), SEM = CFA, raw = FALSE)-model_out$fit
  testthat::expect_equal(abs(fitDiff[1,1]) < 1e-4, TRUE)
  
  fitDiff <- model_out$fit + .5*lambda_*sum(abs(model_out$coefficients[, model_out$pars_pen])) - model_out$optim_fit
  testthat::expect_equal(abs(fitDiff) < 1e-4, TRUE)
  
  
  regularized <- names(pars)[model_out$pars_pen]
  
  aLOOCV2 <- aCV4regsem(lavaanModel = outt,  
                        regsemModel = model_out,
                        k = nrow(HS), 
                        penalty = "lasso",
                        lambda = lambda_)
  
  
  sum(aLOOCV2$leaveOutFits)
  
  exactLOOCV2 <- rep(NA, nrow(HS))
  converged <- exactLOOCV2
  for(i in 1:nrow(HS)){
    outt = cfa(mod, HS[-i,], meanstructure = TRUE)
    model_out = regsem(outt,
                       type="lasso", 
                       pars_pen=c(1:2,6:8),
                       lambda = lambda_, 
                       gradFun = "ram",
                       round = 10)
    converged[i] <- model_out$convergence
    pars <- unlist(model_out$coefficients)
    
    # change labels to those of lavaan
    pars <- pars[labels$regsemLabels]
    names(pars) <- labels$lavaanLabels
    
    CFA <- aCV4SEM:::fit(aCV4SEM:::setParameters(CFA, 
                                                 names(pars), 
                                                 values = pars,
                                                 raw = FALSE))
    
    exactLOOCV2[i] <- aCV4SEM:::computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                                      rawData = as.numeric(HS[i,]),
                                                      impliedMeans = CFA$impliedMeans, 
                                                      impliedCovariance = CFA$impliedCovariance)
    
  }
  
  plot(exactLOOCV2[converged==0], exactLOOCV2[converged==0], type = "l",
       xlab = "exact loocv", ylab = "approximated loocv")
  points(exactLOOCV2[converged==0], aLOOCV2$leaveOutFits[converged==0], col = "red")

  ## reset to test LASSO with aCV
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = FALSE)
  
  # increase to > 25
  cv.out = cv_regsem(outt,type="lasso", 
                     fit.ret = "AIC",
                     metric = "AIC",
                     pars_pen=c(1:2,6:8),
                     n.lambda=50,
                     jump=0.001, 
                     round = 10)
  cv <- aCV4cv_regsem(cvregsemModel = cv.out, lavaanModel = outt, k = nrow(HS), penalty = "lasso")
  
})
