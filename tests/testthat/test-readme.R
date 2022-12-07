test_that("testing readme", {
  testthat::skip_on_cran()
  library(lessSEM)
  
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
  
  regsem <- lasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    # in case of lasso and adaptive lasso, we can specify the number of lambda
    # values to use. lessSEM will automatically find lambda_max and fit
    # models for nLambda values between 0 and lambda_max. For the other
    # penalty functions, lambdas must be specified explicitly
    nLambdas = 50)
  
  # use the plot-function to plot the regularized parameters:
  plot(regsem)
  
  # elements of regsem can be accessed with the @ operator:
  regsem@parameters[1,]
  
  # AIC and BIC:
  AIC(regsem)
  BIC(regsem)
  
  # The best parameters can also be extracted with:
  coef(regsem, criterion = "AIC")
  coef(regsem, criterion = "BIC")
  
  # cross-validation
  cv <- cvLasso(lavaanModel = lavaanModel,
                regularized = paste0("l", 6:15),
                lambdas = seq(0,1,.1),
                standardize = TRUE)
  
  # get best model according to cross-validation:
  coef(cv)
  
  #### Advanced ###
  # Switching the optimizer # 
  # Use the "method" argument to switch the optimizer. The control argument
  # must also be changed to the corresponding function:
  regsemIsta <- lasso(
    lavaanModel = lavaanModel,
    regularized = paste0("l", 6:15),
    nLambdas = 50,
    method = "ista",
    control = controlIsta())
  
  # Note: The results are basically identical:
  regsemIsta@parameters - regsem@parameters
  
})
