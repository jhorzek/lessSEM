# aCV4SEM

aCV4SEM is an R package which provides regularized structural equation modeling (regsem) as well as approximate cross-validation and approximate influence functions building on lavaan. The package is relatively new and you may find more stable implementations of regsem in the R packages regsem and lslx. 

The following features are implemented in aCV4SEM:

- regularization of SEM with ridge, lasso, adaptive lasso, and elastic net (see ?aCV4SEM::regularizeSEM)
- automatic selection of $\lambda$ values for lasso and adaptive lasso (see ?aCV4SEM::regularizeSEM)
- approximate optimization of SEM with custom penalty functions using a BFGS optimizer or Rsolnp (see ?aCV4SEM::regularizeSEMWithCustomPenalty and ?aCV4SEM::regularizeSEMWithCustomPenaltyRsolnp)
- approximate cross-validation (see ?aCV4SEM::aCV4lavaan, ?aCV4SEM::aCV4RegularizedSEM, and ?aCV4SEM::aCV4regularizedSEMWithCustomPenalty)
- approximate influence functions for regularized SEM (see ?aCV4SEM::aI4RegularizedSEM)

# Installation

If you want to install aCV4SEM from GitHub, use the following commands in R:

    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/aCV4SEM")
    

# Example

    library(lavaan) # aCV4SEM builds on lavaan models
    library(aCV4SEM)
    
    ## Approximate leave one out cross-validation for a lavaan model
    ### set up model in lavaan
    HS.model <- ' visual  =~ x1 + x2 + x3
                    textual =~ x4 + x5 + x6
                    speed   =~ x7 + x8 + x9 '
    HS <- HolzingerSwineford1939[,paste0("x",1:9)]
    
    fit <- cfa(HS.model, data = HS, meanstructure = TRUE)
    
    ### approximate cross-validation
    aLOOCV <- aCV4lavaan(lavaanModel = fit, k = nrow(HS))
    
    ### we want to compare this to a true leave one out cross-validation
    exactLOOCV <- rep(NA, nrow(HS))
    for(i in 1:nrow(HS)){
      fit = sem(HS.model, HS[-i,], meanstructure = TRUE)
      exactLOOCV[i] <- aCV4SEM:::computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                                       rawData = as.numeric(HS[i,]),
                                                       impliedMeans = fit@implied$mean[[1]], 
                                                       impliedCovariance = fit@implied$cov[[1]])
    }
    
    # The plot shows the relation between exact and approximate cross-validation.
    # If the points are on the line, the approximate and exact cross-validation
    # produce relatively similar results
    plot(exactLOOCV, exactLOOCV, type = "l",
         xlab = "exact loocv", ylab = "approximated loocv")
    points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
    
    ## Example for regularized model SEM
    
    ### Simulate a data set
    set.seed(123)
    library(mvtnorm)
    # library(semPlot)
    N <- 50
    f <- matrix(rnorm(N, 0, 1), ncol = 1)
    L <- matrix(c(rep(1,5), rep(.2,5), rep(0,5)), nrow = 1)
    y <- matrix(NA, nrow = N, ncol = ncol(L))
    
    covs <- c(rep(1.2-1^2,5), rep(1.2-.8^2,5), rep(1.2-0^2,5))
    
    for(i in 1:N){
      y[i,] <- L*f[i,] +  mvtnorm::rmvnorm(1, sigma = diag(covs))
    }
    
    yNames <- paste0("y", 1:ncol(y))
    colnames(y) <- yNames
    
    ### Fit model with lavaan
    modelSyntax <- paste0('f =~ 1*', yNames[1], ' + ', paste0(yNames[2:length(yNames)], collapse = " + "))
    modelFit = cfa(modelSyntax, y, meanstructure = TRUE)
    
    ### optional: plot model
    # semPlot::semPaths(modelFit)
    
    ### Use lasso regularization for a subset of the loadings.
    regularizedParameterLabels <- paste0("f=~y", 6:15)
    
    regularizedModel <- regularizeSEM(lavaanModel = modelFit,
                                      penalt = "lasso",
                                      regularizedParameterLabels = regularizedParameterLabels,
                                      lambda = seq(0,1,.1)
    )
    
    plot(regularizedModel)
    
    ### approximate leave one out cross-validation:
    aCV <- aCV4regularizedSEM(regularizedModel, k = N)
    plot(aCV)