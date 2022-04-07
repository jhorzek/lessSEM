# aCV4SEM

aCV4SEM provides approximate cross-validation for Structural Equation Models building on lavaan.

# Installation

If you want to install aCV4SEM from GitHub, use the following commands in R:

    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/aCV4SEM")
    

# Example

The following example is adapted from the documentation of ?lavaan::cfa.

    library(lavaan)
    library(aCV4SEM)
    
    HS.model <- ' visual  =~ x1 + x2 + x3
                  textual =~ x4 + x5 + x6
                  speed   =~ x7 + x8 + x9 '
    HS <- HolzingerSwineford1939[,paste0("x",1:9)]
    
    fit <- cfa(HS.model, data = HS, meanstructure = TRUE)
    
    aLOOCV <- approximateCrossValidation(lavaanModel = fit, k = nrow(HS))
    
    
    exactLOOCV <- rep(NA, nrow(HS))
    for(i in 1:nrow(HS)){
      fit = sem(HS.model, HS[-i,], meanstructure = TRUE)
      exactLOOCV[i] <- computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                             rawData = as.numeric(HS[i,]),
                                             impliedMeans = fit@implied$mean[[1]], 
                                             impliedCovariance = fit@implied$cov[[1]])
    }
    
    
    plot(exactLOOCV, exactLOOCV, type = "l",
         xlab = "exact loocv", ylab = "approximated loocv")
    points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
    
    ## Example for regularized model adapted from ?regsem::cv_regsem
    library(regsem)
    # put variables on same scale for regsem
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
    
    aCV <- approximateCrossValidation(lavaanModel = outt, # lavaan model is always required as basis
                                      SEM = cv.out, 
                                      k = nrow(HS), # leave one out cross-validation
                                      penalty = "lasso",
                                      eps = 1e-3)
    cvFits <- apply(aCV$approximateCV,2,sum)
    par(mfrow = c(1,2))
    plot(x = aCV$lambda, 
         y = cvFits,
         xlab = "lambda", 
         ylab = "aCV", 
         type = "l")
    abline(v = aCV$lambda[which(cvFits == min(cvFits))])
    
    plot(x = aCV$lambda, 
         y = cv.out$fits[,"AIC"],
         xlab = "lambda", 
         ylab = "AIC", 
         type = "l")
    abline(v = aCV$lambda[which(cv.out$fits[,"AIC"] == min(cv.out$fits[,"AIC"]))])
    par(mfrow = c(1,1))
