# aCV4SEM

aCV4SEM provides approximate cross-validation for Structural Equation Models building no lavaan.

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
