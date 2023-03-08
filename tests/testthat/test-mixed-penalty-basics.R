test_that("testing mixedPenalty-basic", {
  testthat::skip_on_cran()
  library(lessSEM)
  set.seed(123)
  N <- 100
  l1 <- 1; l2 <- .2; l3 <- 0;
  v1 <- .2; v2 <- .8; v3 <- 1
  
  f <- matrix(stats::rnorm(N, 0, 1), ncol = 1)
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
  
  # fit regularizedSEM
  lambdas <- seq(0,1,length.out=5)
  regularized = paste0("f=~y",6:ncol(y))
  rsemIsta <- lasso(lavaanModel = modelFit, 
                    regularized = regularized, 
                    lambdas = lambdas,
                    method = "ista",
                    control = controlIsta()
  )
  
  # fit with mixedPenalty
  mixed <- mixedPenalty(modelFit) |>
    addLasso(regularized = regularized,
             lambdas = lambdas) |>
    fit()
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] - mixed@parameters[,rsemIsta@regularized]) < .002), TRUE)
  coef(mixed, criterion = "AIC")
  coef(mixed, criterion = "BIC")
  
  # fit regularizedSEM
  rsemIsta <- cappedL1(lavaanModel = modelFit, 
                       regularized = regularized, 
                       lambdas = lambdas,
                       thetas = .2
  )
  
  # fit with mixedPenalty
  mixed <- mixedPenalty(modelFit) |>
    addCappedL1(regularized = regularized,
                lambdas = lambdas,
                thetas = .2) |>
    fit()
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] - mixed@parameters[,rsemIsta@regularized]) < .002), TRUE)
  coef(mixed, criterion = "AIC")
  coef(mixed, criterion = "BIC")
  
  # fit regularizedSEM
  rsemIsta <- lsp(lavaanModel = modelFit, 
                  regularized = regularized, 
                  lambdas = lambdas,
                  thetas = .2
  )
  
  # fit with mixedPenalty
  mixed <- mixedPenalty(modelFit) |>
    addLsp(regularized = regularized,
           lambdas = lambdas,
           thetas = .2) |>
    fit()
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] - mixed@parameters[,rsemIsta@regularized]) < .002), TRUE)
  coef(mixed, criterion = "AIC")
  coef(mixed, criterion = "BIC")
  
  # fit regularizedSEM
  rsemIsta <- mcp(lavaanModel = modelFit, 
                  regularized = regularized, 
                  lambdas = lambdas,
                  thetas = .2
  )
  
  # fit with mixedPenalty
  mixed <- mixedPenalty(modelFit) |>
    addMcp(regularized = regularized,
           lambdas = lambdas,
           thetas = .2) |>
    fit()
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] - mixed@parameters[,rsemIsta@regularized]) < .002), TRUE)
  coef(mixed, criterion = "AIC")
  coef(mixed, criterion = "BIC")
  
  # fit regularizedSEM
  rsemIsta <- scad(lavaanModel = modelFit, 
                   regularized = regularized, 
                   lambdas = lambdas,
                   thetas = 2.2
  )
  
  # fit with mixedPenalty
  mixed <- mixedPenalty(modelFit) |>
    addScad(regularized = regularized,
            lambdas = lambdas,
            thetas = 2.2) |>
    fit()
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] - mixed@parameters[,rsemIsta@regularized]) < .002), TRUE)
  coef(mixed, criterion = "AIC")
  coef(mixed, criterion = "BIC")
  
  set.seed(123)
  # Identical to regsem, lessSEM builds on the lavaan
  # package for model specification. The first step
  # therefore is to implement the model in lavaan.
  
  dataset <- simulateExampleData()
  
  lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"
  
  lavaanModel <- lavaan::sem(lavaanSyntax,
                             data = dataset,
                             meanstructure = TRUE,
                             std.lv = TRUE)
  
  # Optional: Plot the model
  # semPlot::semPaths(lavaanModel, 
  #                   what = "est",
  #                   fade = FALSE)
  
  # We can add multiple penalties as follows:
  mixedPenaltyGlmnet <- lavaanModel |>
    # create template for regularized model with mixed penalty:
    mixedPenalty(method = "glmnet", control = controlGlmnet()) |>
    # add lasso penalty on loadings l6 - l10:
    addElasticNet(regularized = paste0("l", 6:10), 
                  lambdas = seq(0,1,length.out=3),
                  alphas = 1) |>
    # add penalty on loadings l11 - l15:
    addElasticNet(regularized = paste0("l", 11:15), 
                  lambdas = seq(0,1,length.out=3),
                  alphas = 1) |>
    fit()
  
  
  # We can add multiple penalties as follows:
  mixedPenaltyIsta <- lavaanModel |>
    # create template for regularized model with mixed penalty:
    mixedPenalty() |>
    # add lasso penalty on loadings l6 - l10:
    addLasso(regularized = paste0("l", 6:10), 
             lambdas = seq(0,1,length.out=3)) |>
    # add penalty on loadings l11 - l15:
    addLasso(regularized = paste0("l", 11:15), 
             lambdas = seq(0,1,length.out=3)) |>
    fit()
  
  testthat::expect_equal(all(abs(coef(mixedPenaltyGlmnet)@estimates - coef(mixedPenaltyIsta)@estimates) < .001), TRUE)
  
  # We can add multiple penalties as follows:
  mixedPenaltyGlmnet <- lavaanModel |>
    # create template for regularized model with mixed penalty:
    mixedPenalty(method = "glmnet", control = controlGlmnet()) |>
    # add lasso penalty on loadings l6 - l10:
    addLasso(regularized = paste0("l", 6:10), 
                  lambdas = seq(0,1,length.out=3)) |>
    # add penalty on loadings l11 - l15:
    addLasso(regularized = paste0("l", 11:15), 
                  lambdas = seq(0,1,length.out=3)) |>
    fit()
  
  testthat::expect_equal(all(abs(coef(mixedPenaltyGlmnet)@estimates - coef(mixedPenaltyIsta)@estimates) < .001), TRUE)
  
  testthat::expect_error(
    # We can add multiple penalties as follows:
    mixedPenaltyGlmnet <- lavaanModel |>
      # create template for regularized model with mixed penalty:
      mixedPenalty(method = "glmnet", control = controlGlmnet()) |>
      # add lasso penalty on loadings l6 - l10:
      addElasticNet(regularized = paste0("l", 6:10), 
                    lambdas = seq(0,1,length.out=3),
                    alphas = 1) |>
      # add penalty on loadings l11 - l15:
      addScad(regularized = paste0("l", 11:15), 
               lambdas = seq(0,1,length.out=3), 
              thetas = 2.4) |>
      fit()
    
  )
  
  testthat::expect_error(
    # We can add multiple penalties as follows:
    mixedPenaltyIsta <- lavaanModel |>
      # create template for regularized model with mixed penalty:
      mixedPenalty() |>
      # add lasso penalty on loadings l6 - l10:
      addElasticNet(regularized = paste0("l", 6:10), 
                    lambdas = seq(0,1,length.out=3),
                    alphas = 1) |>
      # add scad penalty on loadings l11 - l15:
      addScad(regularized = paste0("l", 11:15), 
              lambdas = seq(0,1,length.out=3), 
              thetas = 2.4) |>
      fit()
  )
  
})
