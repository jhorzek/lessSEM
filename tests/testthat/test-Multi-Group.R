test_that("testing multi-group SEM", {
  testthat::skip_on_cran()
  library(lessSEM)
  set.seed(123)
  # Shortened example from https://lavaan.ugent.be/tutorial/groups.html
  
  # model 1: all parameters are the same
  HS.model <- ' visual  =~ c(l1,l1)*x1 + c(l2,l2)*x2 + c(l3,l3)*x3
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
  
  groups <- c("Pasteur", "Grant-White")
  lavaanModels <- c()
  for(i in 1:length(groups)){
    HS.model_i <- paste0(' 
    visual  =~ l1*x1 + l2*x2 + l3*x3
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
  
  regsem <- bfgs(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels)
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  # same model, but now with group-specific loadings:
  
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
  
  regsem <- bfgs(
    # pass the fitted lavaan models
    lavaanModel = lavaanModels)
  
  testthat::expect_equal(any(abs(unlist(regsem@parameters[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
  
  # same thing, but using a transformation:
  
  transformation <- "
parameters: l1_1, l1_2, delta1, l2_1, l2_2, delta2, l3_1, l3_2, delta3
l1_2 = l1_1 + delta1;
l2_2 = l2_1 + delta2;
l3_2 = l3_1 + delta3;
"
  regsemTransform <- bfgs(
    # pass the fitted lavaan model
    lavaanModel = lavaanModels,
    modifyModel = modifyModel(transformations = transformation))
  
  testthat::expect_equal(any(abs(unlist(cbind(regsemTransform@parameters, regsemTransform@transformations)[1,names(getLavaanParameters(fit))]) - getLavaanParameters(fit)) > .001),
                         FALSE)
})
