test_that("testing regularized multi-group SEM", {
  testthat::skip_on_cran()
  library(lessSEM)
  set.seed(123)
  # Shortened example from https://lavaan.ugent.be/tutorial/groups.html
  
  groups <- c("Pasteur", "Grant-White")
  # group-specific loadings:
  
  HS.model <- ' visual  =~ c(l1_1,l1_2)*x1 + c(l2_1,l2_2)*x2 + c(l3_1,l3_2)*x3
x1 ~~ c(v1,v1)*x1
x2 ~~ c(v2,v2)*x2
x3 ~~ c(v3,v3)*x3

visual ~~ c(lv1,lv1)*visual
x1 ~ c(m1,m1)*1
x2 ~ c(m2,m2)*1
x3 ~ c(m3,m3)*1
'
  
  fit <- cfa(HS.model, 
             data = HolzingerSwineford1939, 
             group = "school",
             std.lv = TRUE)
  
  lavaanModels <- c()
  for(i in 1:length(groups)){
    HS.model_i <- paste0(' 
    visual  =~ l1_',i, '*x1 + l2_',i, '*x2 + l3_',i, '*x3
    x1 ~~ v1*x1
    x2 ~~ v2*x2
    x3 ~~ v3*x3
    
    visual ~~ lv1*visual
    x1 ~ m1*1
    x2 ~ m2*1
    x3 ~ m3*1')
    lavaanModels <- c(lavaanModels,
                      lavaan::sem(HS.model_i,
                                  data = HolzingerSwineford1939[HolzingerSwineford1939$school == groups[i],],
                                  meanstructure = TRUE,
                                  std.lv = TRUE)
    )
  }
  
  ### lasso ####
  regsem <- lasso(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = "l1_1",
    nLambdas = 5)
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[nrow(regsem@parameters),names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### adaptive lasso ####
  regsem <- adaptiveLasso(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = "l1_1",
    nLambdas = 5)
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[nrow(regsem@parameters),names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### elastic net ####
  regsem <- elasticNet(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = "l1_1",
    lambdas = seq(0,1,length.out = 3),
    alphas = seq(0,1,.25))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### scad ####
  regsem <- scad(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = "l1_1",
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(2.5,3,.5))
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### mcp ####
  regsem <- mcp(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = "l1_1",
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(.5,1,.5))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### lsp ####
  regsem <- lsp(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = "l1_1",
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(.5,1,5))
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### cappedL1 ####
  regsem <- cappedL1(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = "l1_1",
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(.5,1,.5))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  # same thing, but using a transformation:
  
  transformation <- "
parameters: l1_1, l1_2, delta1, l2_1, l2_2, delta2, l3_1, l3_2, delta3
l1_2 = l1_1 + delta1;
l2_2 = l2_1 + delta2;
l3_2 = l3_1 + delta3;
"
  
  ### lasso ####
  
  regsem <- lasso(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = c("delta1", "delta2", "delta3"),
    nLambdas = 5,
    modifyModel = modifyModel(transformations = transformation))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(cbind(regsem@parameters, regsem@transformations)[nrow(regsem@parameters),names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### adaptive lasso ####
  regsem <- adaptiveLasso(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = c("delta1", "delta2", "delta3"),
    nLambdas = 5,
    modifyModel = modifyModel(transformations = transformation))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(cbind(regsem@parameters, regsem@transformations)[nrow(regsem@parameters),names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### elastic net ####
  regsem <- elasticNet(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = c("delta1", "delta2", "delta3"),
    lambdas = seq(0,1,length.out = 3),
    alphas = seq(0,1,.5),
    modifyModel = modifyModel(transformations = transformation))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(cbind(regsem@parameters, regsem@transformations)[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### scad ####
  regsem <- scad(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = c("delta1", "delta2", "delta3"),
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(2.5,3,.5),
    modifyModel = modifyModel(transformations = transformation))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(cbind(regsem@parameters, regsem@transformations)[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### mcp ####
  regsem <- mcp(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = c("delta1", "delta2", "delta3"),
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(.5,1,.5),
    modifyModel = modifyModel(transformations = transformation))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(cbind(regsem@parameters, regsem@transformations)[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### lsp ####
  regsem <- lsp(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = c("delta1", "delta2", "delta3"),
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(.5,1,.5),
    modifyModel = modifyModel(transformations = transformation))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(cbind(regsem@parameters, regsem@transformations)[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  ### cappedL1 ####
  regsem <- cappedL1(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels,
    regularized = c("delta1", "delta2", "delta3"),
    lambdas = seq(0,1,length.out = 3),
    thetas = seq(.5,1,.5),
    modifyModel = modifyModel(transformations = transformation))
  
  AIC(regsem)
  BIC(regsem)
  
  testthat::expect_equal(any(abs(unlist(cbind(regsem@parameters, regsem@transformations)[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
})
