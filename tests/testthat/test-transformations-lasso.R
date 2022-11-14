test_that("testing elasticNet-lasso-with-transformation", {
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  model <- ' 
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
  
  modelFit <- sem(model, 
               data = PoliticalDemocracy, 
               meanstructure = TRUE, 
               do.fit = TRUE)
  
  # let's define some transformations to test for measurement invariance:
  transformations <- "
  parameters: a,b,c,d,e,f,deltaA, deltaB, deltaC
  d = a + deltaA;
  e = b + deltaB;
  f = c + deltaC;
  "
  
  lambdas <- seq(0,.4,length.out = 10)
  
  rsemIsta <- lasso(lavaanModel = modelFit, 
                    regularized = c("deltaA", "deltaB", "deltaC"), 
                    lambdas = lambdas,
                    modifyModel = modifyModel(transformations = transformations)
  )
  
  plot(rsemIsta)
  coef(rsemIsta)
  coef(rsemIsta, criterion = "AIC")
  coef(rsemIsta, criterion = "BIC")
  
  rsemGlmnet <- lasso(lavaanModel = modelFit, 
                      regularized = c("deltaA", "deltaB", "deltaC"),
                      lambdas = lambdas,
                      method = "glmnet",
                      control = controlGlmnet(),
                      modifyModel = modifyModel(transformations = transformations)
  )
  plot(rsemGlmnet)
  coef(rsemGlmnet)
  coef(rsemGlmnet, criterion = "AIC")
  coef(rsemGlmnet, criterion = "BIC")
  
  testthat::expect_equal(all(abs(coef(rsemGlmnet) - coef(rsemIsta)) < 1e-2), TRUE)
  
  # test multi-line transformations
  transformations <- "
  parameters: a,b,c,d,e,f,
  deltaA, deltaB, 
  deltaC
  d = a + 
  deltaA;
  e = b + deltaB;
  f = c + ((
  deltaC
  )
  );
  "
  rsemGlmnet2 <- lasso(lavaanModel = modelFit, 
                      regularized = c("deltaA", "deltaB", "deltaC"),
                      lambdas = lambdas,
                      method = "glmnet",
                      control = controlGlmnet(),
                      modifyModel = modifyModel(transformations = transformations)
  )
  testthat::expect_equal(all(abs(coef(rsemGlmnet) - coef(rsemGlmnet2)) < 1e-4), TRUE)
  
  # test multi-line transformations and starting values
  transformations <- "
  parameters: a,b,c,d,e,f,
  deltaA, deltaB, 
  deltaC
  start: a = 1, b = 2, 
  c=3, 
  deltaB = 1.1, 
    deltaC = .2
  d = a + deltaA;
  e = b + deltaB;
  f = c + deltaC;
  "
  rsemGlmnet3 <- lasso(lavaanModel = modelFit, 
                       regularized = c("deltaA", "deltaB", "deltaC"),
                       lambdas = lambdas,
                       method = "glmnet",
                       control = controlGlmnet(),
                       modifyModel = modifyModel(transformations = transformations)
  )
  testthat::expect_equal(all(abs(coef(rsemGlmnet) - coef(rsemGlmnet3)) < 1e-2), TRUE)
  
  # set automatic lambda:
  transformations <- "
  parameters: a,b,c,d,e,f,deltaA, deltaB, deltaC
  d = a + deltaA;
  e = b + deltaB;
  f = c + deltaC;
  "
  rsem2 <- lessSEM::lasso(lavaanModel = modelFit, 
                          regularized = c("deltaA", "deltaB", "deltaC"),
                          nLambdas = 10,
                          modifyModel = modifyModel(transformations = transformations))
  testthat::expect_equal(all(apply(rsem2@parameters[,c("deltaA", "deltaB", "deltaC")] == 0,2,sum) > 0), TRUE)
  
  rsem2 <- lessSEM::lasso(lavaanModel = modelFit, 
                          regularized = c("deltaA", "deltaB", "deltaC"),
                          nLambdas = 10,
                          modifyModel = modifyModel(transformations = transformations),
                          reverse = FALSE)
  testthat::expect_equal(all(apply(rsem2@parameters[,c("deltaA", "deltaB", "deltaC")] == 0,2,sum) > 0), TRUE)
  
  # cross-validation
  rsemGlmnet <- cvLasso(lavaanModel = modelFit, 
                      regularized = c("deltaA", "deltaB", "deltaC"),
                      lambdas = lambdas,
                      method = "glmnet",
                      control = controlGlmnet(),
                      modifyModel = modifyModel(transformations = transformations),
                      k = 5
  )
  rsemGlmnet@parameters
  
  # adaptive lasso
  rsemGlmnet <- adaptiveLasso(lavaanModel = modelFit, 
                        regularized = c("deltaA", "deltaB", "deltaC"),
                        lambdas = lambdas,
                        method = "glmnet",
                        control = controlGlmnet(),
                        modifyModel = modifyModel(transformations = transformations)
  )
  
  # compare weights to mles:
  MLEs <- bfgs(lavaanModel = modelFit, 
                              modifyModel = modifyModel(transformations = transformations)
  )
  testthat::expect_equal(all(abs(rsemGlmnet@weights[c("deltaA", "deltaB", "deltaC")]^(-1) - 
                               abs(MLEs@parameters[,c("deltaA", "deltaB", "deltaC")])) < 1e-3),TRUE)
  
  
})
