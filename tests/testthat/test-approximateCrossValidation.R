test_that("approximate cross validation works", {
  
  library(lavaan)
  library(regsem)
  
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = TRUE)
  
  CFA <- SEMFromLavaan(outt, rawData = HS)
  CFA <- fit(CFA)
  
  hessian <- computeHessian(SEM = CFA, fitfunction = minus2LogLikelihood)
  
  aLOOCV <- approximateCrossValidation(SEM = CFA, k = nrow(HS), hessian = hessian, individualFitfunction = individualMinus2LogLikelihood)
  
  
  exactLOOCV <- rep(NA, nrow(HS))
  for(i in 1:nrow(HS)){
    outt = cfa(mod, HS[-i,], meanstructure = TRUE)
    exactLOOCV[i] <- computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                           rawData = as.numeric(HS[i,]),
                                           expectedMeans = outt@implied$mean[[1]], 
                                           expectedCovariance = outt@implied$cov[[1]])
  }
  
  testthat::expect_equal(abs(sum(aLOOCV$leaveOutFits) - sum(exactLOOCV)) < 2, TRUE)
  
  plot(exactLOOCV, exactLOOCV, type = "l",
       xlab = "exact loocv", ylab = "approximated loocv")
  points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
  
  # now try regularized models
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = TRUE)
  
  CFA <- SEMFromLavaan(outt, rawData = HS)
  CFA <- fit(CFA)
  
  lambda_ <- 0.05
  model_out = regsem(outt,
                     type="lasso", 
                     pars_pen=c(1:2,6:8),
                     lambda = lambda_,
                     gradFun = "ram",
                     round = 10)
  CFA <- setParameters(CFA, labels = names(model_out$coefficients), 
                       values = unlist(model_out$coefficients), 
                       labelsFrom = "regsem", 
                       raw = FALSE)
  likelihoodRatioFitRegsem(par = getParameters(CFA), SEM = CFA, raw = FALSE)
  model_out$fit
  
  model_out$fit + .5*lambda_*sum(abs(model_out$coefficients[, model_out$pars_pen]))
  model_out$optim_fit
  
  regularized <- names(getParameters(CFA))[model_out$pars_pen]
  
  penaltyFunction <- function(par, lambda, regularizedParameters, Npen){
    return(lambda*sum(sqrt((par[regularizedParameters])^2 + 1e-4)))
  }
  
  penaltyScoreFunction <- function(par, lambda, regularized, Npen){
    scores <- matrix(0, nrow = Npen, ncol = length(par))
    colnames(scores) <- names(par)
    for(p in regularized){
      scores[,p] <- (1/Npen)*lambda*par[p]/(sqrt(par[p]^2 + 1e-4))
    }
    return(scores)
  }
  
  Npen <- nrow(HS)
  
  
  aLOOCV2 <- approximateCrossValidation(SEM = CFA, 
                                        k = Npen, 
                                        penaltyScoreFunction = penaltyScoreFunction, 
                                        penaltyFunction = penaltyFunction, 
                                        lambda = Npen*lambda_, 
                                        regularized = regularized, 
                                        Npen = Npen)

  
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
    CFA <- fit(setParameters(CFA, 
                             names(model_out$coefficients), 
                             values = unlist(model_out$coefficients),
                             labelsFrom = "regsem", 
                             raw = FALSE))
    
    exactLOOCV2[i] <- computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                           rawData = as.numeric(HS[i,]),
                                           expectedMeans = CFA$model$expected$means, 
                                           expectedCovariance = CFA$model$expected$covariance)
    
  }
  
  plot(exactLOOCV2[converged==0], exactLOOCV2[converged==0], type = "l",
       xlab = "exact loocv", ylab = "approximated loocv")
  points(exactLOOCV2[converged==0], aLOOCV2$leaveOutFits[converged==0], col = "red")
  #testthat::expect_equal(abs(sum(aLOOCV2$leaveOutFits) - sum(exactLOOCV)) < 2, TRUE)
  
})
