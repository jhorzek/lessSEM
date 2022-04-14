test_that("approximate cross validation works", {
  
  library(lavaan)
  library(regsem)
  
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = TRUE)
  
  aLOOCV <- aCV4lavaan(lavaanModel = outt, k = nrow(HS), raw = FALSE)
  
  exactLOOCV <- rep(NA, nrow(HS))
  for(i in 1:nrow(HS)){
    outt = cfa(mod, HS[-i,], meanstructure = TRUE)
    exactLOOCV[i] <- computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                           rawData = as.numeric(HS[i,]),
                                           impliedMeans = outt@implied$mean[[1]], 
                                           impliedCovariance = outt@implied$cov[[1]])
  }
  
  testthat::expect_equal(max(abs(aLOOCV$leaveOutFits -exactLOOCV)) < .1, TRUE)
  
  plot(exactLOOCV, exactLOOCV, type = "l",
       xlab = "exact loocv", ylab = "approximated loocv")
  points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
  
  # now try regularized models
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = TRUE)
  
  CFA <- SEMFromLavaan(lavaanModel = outt)
  CFA <- fit(CFA)
  
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
  
  CFA <- setParameters(CFA, 
                       labels = names(pars), 
                       values = pars, 
                       raw = FALSE)
  .5*(1/nrow(HS))*likelihoodRatioFit(par = getParameters(CFA), SEM = CFA, raw = FALSE)
  model_out$fit
  
  model_out$fit + .5*lambda_*sum(abs(model_out$coefficients[, model_out$pars_pen]))
  model_out$optim_fit
  
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
    
    CFA <- fit(setParameters(CFA, 
                             names(pars), 
                             values = pars,
                             raw = FALSE))
    
    exactLOOCV2[i] <- computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                            rawData = as.numeric(HS[i,]),
                                            impliedMeans = CFA$impliedMeans, 
                                            impliedCovariance = CFA$impliedCovariance)
    
  }
  
  plot(exactLOOCV2[converged==0], exactLOOCV2[converged==0], type = "l",
       xlab = "exact loocv", ylab = "approximated loocv")
  points(exactLOOCV2[converged==0], aLOOCV2$leaveOutFits[converged==0], col = "red")
  #testthat::expect_equal(abs(sum(aLOOCV2$leaveOutFits) - sum(exactLOOCV)) < 2, TRUE)
  
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
  cvregsemPars <- cvregsem2LavaanParameters(cvregsemModel = cv.out, lavaanModel = outt)
  
  regularizedParameterLabels <- paste0("f=~x", c(2:3, 7:9))
  optLasso2 <- regularizeSEM(lavaanModel = outt,
                             regularizedParameterLabels = regularizedParameterLabels, 
                             lambda = cv.out$fits[,"lambda"], 
                             penalty = "lasso")
  testthat::expect_equal(all(round(wideResults(optLasso2)[,colnames(cvregsemPars)] - cvregsemPars,1) == 0), TRUE)
  matplot(optLasso2@inputArguments$lambdas, wideResults(optLasso2), type = "l")
  matplot(optLasso2@inputArguments$lambdas, cvregsemPars, type = "l")
  cv <- aCV4regularizedSEM(regularizedSEM = optLasso2, k = nrow(HS))
  coef(cv)
  plot(cv)
  
  ## Test ridge
  optLasso3 <- regularizeSEM(lavaanModel = outt,
                             regularizedParameterLabels = regularizedParameterLabels, 
                             lambda = cv.out$fits[,"lambda"], 
                             penalty = "ridge")
  plot(optLasso3)
  cv <- aCV4regularizedSEM(regularizedSEM = optLasso3, k = nrow(HS))
  coef(cv)
  plot(cv)

  ## Test adaptive lasso
  optLasso4 <- regularizeSEM(lavaanModel = outt,
                             regularizedParameterLabels = regularizedParameterLabels, 
                             lambda = cv.out$fits[,"lambda"], 
                             penalty = "adaptiveLasso")
  plot(optLasso4)
  cv <- aCV4regularizedSEM(regularizedSEM = optLasso4, k = nrow(HS))
  coef(cv)
  plot(cv)

  
  ## Test elastic net
  optLasso5 <- regularizeSEM(lavaanModel = outt,
                             regularizedParameterLabels = regularizedParameterLabels, 
                             lambda = c(0,.2,.5), 
                             alpha = c(0,.1),
                             penalty = "elasticNet")
  cv <- aCV4regularizedSEM(regularizedSEM = optLasso5, k = nrow(HS))
  coef(cv)
  plot(cv, alpha = .1)
})
