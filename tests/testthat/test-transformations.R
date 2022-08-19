test_that("testing transformations", {
  library(lavaan)
  library(lessSEM)
  model1 <- ' 
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
  
  model <- sem(model1, 
               data = PoliticalDemocracy, 
               meanstructure = TRUE, 
               do.fit = FALSE)
  
  # let's test a simple equality constraint:
  transformations <- "
  parameters: a,b,c,d,e,f
  d = a
  e = b
  f = c
  "
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model, 
                                  whichPars = "start",
                                  fit = TRUE, 
                                  addMeans = TRUE,
                                  transformations = transformations)
  SEM$fit()
  
  parameters <- lessSEM:::.getParameters(SEM, raw = TRUE)
  gradients <- lessSEM:::.getGradients(SEM, raw = TRUE)
  
  # these gradients should be identical to those of a model where the
  # equality constraint was present from the get-go
  model2 <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

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
  model2 <- sem(model2, 
                data = PoliticalDemocracy, 
                meanstructure = TRUE,
                do.fit = FALSE
  )
  
  SEMNoTransformations <- lessSEM:::.SEMFromLavaan(lavaanModel = model2, 
                                                   whichPars = "start",
                                                   fit = TRUE, 
                                                   addMeans = TRUE)
  SEMNoTransformations <- lessSEM:::.setParameters(SEMNoTransformations, 
                                                   names(parameters),
                                                   parameters, 
                                                   TRUE)
  
  testthat::expect_equal(abs(SEM$fit()- SEMNoTransformations$fit()) < 1e-8,TRUE)
  
  gradientsNoTransformations <- lessSEM:::.getGradients(SEMNoTransformations, raw = TRUE)
    
  testthat::expect_equal(all(abs(gradients - gradientsNoTransformations[names(gradients)]) < 1e-6), TRUE)
  
  # Now, let's test a more complicated transformation:
    
    transformations <- "
  parameters: a,b,c,d,e,f,deltaA, deltaB, deltaC
  d = a + deltaA
  e = b + deltaB
  f = c + deltaC
  "
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model, 
                                  whichPars = "start",
                                  fit = FALSE, 
                                  addMeans = TRUE,
                                  transformations = transformations)
  
  params <- lessSEM:::.getParameters(SEM, raw = TRUE, transformations = TRUE)
  testthat::expect_equal((params["d"] - params["a"] - params["deltaA"]) == 0,c("d" = TRUE))
  testthat::expect_equal((params["e"] - params["b"] - params["deltaB"]) == 0,c("e" = TRUE))
  testthat::expect_equal((params["f"] - params["c"] - params["deltaC"]) == 0,c("f" = TRUE))
  
  params <- lessSEM:::.getParameters(SEM, raw = TRUE, transformations = FALSE)
  params["deltaA"] <- .2
  params["deltaB"] <- -.1
  params["deltaC"] <- .3
  SEM <- lessSEM:::.setParameters(SEM, labels = names(params), values = params, raw = TRUE)
  params2 <- lessSEM:::.getParameters(SEM, raw = TRUE, transformations = TRUE)
  
  testthat::expect_equal(abs(params2["d"] - params2["a"] - params2["deltaA"])< 1e-10,c("d" = TRUE))
  testthat::expect_equal(abs(params2["e"] - params2["b"] - params2["deltaB"])< 1e-10,c("e" = TRUE))
  testthat::expect_equal(abs(params2["f"] - params2["c"] - params2["deltaC"])< 1e-10,c("f" = TRUE))
  
})
