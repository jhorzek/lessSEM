test_that("testing approximate lasso", {
  library(lavaan)
  library(linr)
  set.seed(123)
  N <- 100
  l1 <- 1; l2 <- .2; l3 <- 0;
  v1 <- .2; v2 <- .8; v3 <- 1
  
  f <- matrix(rnorm(N, 0, 1), ncol = 1)
  L <- matrix(c(rep(l1,5), rep(l2,5), rep(l3,15)), nrow = 1)
  y <- matrix(NA, nrow = N, ncol = ncol(L))
  
  covs <- c(rep(v1,5), rep(v2,5), rep(v3,15))
  
  for(i in 1:N){
    y[i,] <- L*f[i,] +  mvtnorm::rmvnorm(1, sigma = diag(covs))
  }
  
  yNames <- paste0("y", 1:ncol(y))
  colnames(y) <- yNames
  
  modelSyntax <- paste0('f =~ 1*', yNames[1], ' + ', paste0(yNames[2:length(yNames)], collapse = " + "))
  modelFit = cfa(modelSyntax, y, meanstructure = TRUE)
  
  
  # regularize 
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  rsem <- regularizeSEM(lavaanModel = modelFit, 
                        regularizedParameterLabels = regularizedLavaan,
                        penalty = "lasso", 
                        nLambdas = 30)
  
  lambdas <- rsem@fits$lambda 
  
  penaltyFunctionArguments <- list(
    regularizedParameterLabels = regularizedLavaan,
    eps = 1e-7 # note: the optimization works better, if eps is NOT set to 0! This is especially true for aLOOCV
  )
  tuningParameters <- data.frame("lambda" = lambdas)
  
  apprRegsem <- regularizeSEMWithCustomPenalty(lavaanModel = modelFit, 
                                               individualPenaltyFunction = smoothLASSO, 
                                               individualPenaltyFunctionGradient = smoothLASSOGradient, 
                                               individualPenaltyFunctionHessian = smoothLASSOHessian, 
                                               tuningParameters = tuningParameters, 
                                               penaltyFunctionArguments = penaltyFunctionArguments)
  
  parameterDifference <- apprRegsem@parameters[,rsem@parameterLabels] - rsem@parameters[,rsem@parameterLabels]
  matplot(abs(parameterDifference[,regularizedLavaan]), type  ="l")
  testthat::expect_equal(max(abs(parameterDifference[,regularizedLavaan])) < .03, TRUE)
  
  ## with approximate Hessian
  apprRegsem2 <- regularizeSEMWithCustomPenalty(lavaanModel = modelFit, 
                                                individualPenaltyFunction = smoothLASSO, 
                                                individualPenaltyFunctionGradient = smoothLASSOGradient,
                                                tuningParameters = tuningParameters, 
                                                penaltyFunctionArguments = penaltyFunctionArguments)
  
  parameterDifference2 <- apprRegsem2@parameters[,rsem@parameterLabels] - rsem@parameters[,rsem@parameterLabels]
  matplot(abs(parameterDifference2[,regularizedLavaan]), type  ="l")
  testthat::expect_equal(max(abs(parameterDifference2[,regularizedLavaan])) < .03, TRUE)
  
  # approximate gradient AND Hessian
  apprRegsem3 <- regularizeSEMWithCustomPenalty(lavaanModel = modelFit, 
                                                individualPenaltyFunction = smoothLASSO,
                                                tuningParameters = tuningParameters, 
                                                penaltyFunctionArguments = penaltyFunctionArguments)
  
  parameterDifference3 <- apprRegsem3@parameters[,rsem@parameterLabels] - rsem@parameters[,rsem@parameterLabels]
  matplot(abs(parameterDifference3[,regularizedLavaan]), type  ="l")
  testthat::expect_equal(max(abs(parameterDifference3[,regularizedLavaan])) < .03, TRUE)
  
  # with Rsolnp:
  apprRegsem4 <- regularizeSEMWithCustomPenaltyRsolnp(lavaanModel = modelFit, 
                                                      individualPenaltyFunction = smoothLASSO, 
                                                      tuningParameters = tuningParameters, 
                                                      penaltyFunctionArguments = penaltyFunctionArguments)
  parameterDifference4 <- apprRegsem4@parameters[,rsem@parameterLabels] - rsem@parameters[,rsem@parameterLabels]
  matplot(abs(parameterDifference4[,regularizedLavaan]), type  ="l")
  testthat::expect_equal(max(abs(parameterDifference4[,regularizedLavaan])) < .03, TRUE)
  
  # test approximate cross-validation
  aCV1 <- aCV4regularizedSEMWithCustomPenalty(regularizedSEMWithCustomPenalty = apprRegsem, k = N)
  print(aCV1)
  coef(aCV1)
  plot(aCV1)
  
  AICs <- AIC(apprRegsem, penalizedParameterLabels = apprRegsem@inputArguments$penaltyFunctionArguments$regularizedParameterLabels, zeroThreshold = 1e-4)
  plot(AICs$lambda, AICs$AIC, type = "l")
  BICs <- BIC(apprRegsem, penalizedParameterLabels = apprRegsem@inputArguments$penaltyFunctionArguments$regularizedParameterLabels, zeroThreshold = 1e-4)
  plot(BICs$lambda, BICs$BIC, type = "l")
})
