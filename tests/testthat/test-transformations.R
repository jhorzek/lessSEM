test_that("testing transformations", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("OpenMx")
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
  d = a;
  e = b;
  f = c;
  "
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model, 
                                  whichPars = "start",
                                  fit = TRUE, 
                                  addMeans = TRUE,
                                  transformations = transformations)
  SEM$fit()
  
  parameters <- lessSEM:::.getParameters(SEM, raw = TRUE)
  gradients <- lessSEM:::.getGradients(SEM, raw = TRUE)
  Hessian <- lessSEM:::.getHessian(SEM, raw = TRUE)
  
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
  parameters_ <- parameters[!names(parameters) %in% c("d", "e", "f")] 
  SEMNoTransformations <- lessSEM:::.setParameters(SEMNoTransformations, 
                                                   names(parameters_),
                                                   parameters_, 
                                                   TRUE)
  
  testthat::expect_equal(abs(SEM$fit()- SEMNoTransformations$fit()) < 1e-8,TRUE)
  
  gradientsNoTransformations <- lessSEM:::.getGradients(SEMNoTransformations, raw = TRUE)
  HessianNoTransformations <- lessSEM:::.getHessian(SEM, raw = TRUE)
  
  testthat::expect_equal(all(abs(gradients[names(gradientsNoTransformations)] -
                                   gradientsNoTransformations[names(gradientsNoTransformations)]) < 1e-6), TRUE)
  
  testthat::expect_equal(all(abs(Hessian[names(gradientsNoTransformations), names(gradientsNoTransformations)] - 
                                   HessianNoTransformations[names(gradientsNoTransformations), names(gradientsNoTransformations)]) < 1e-6), TRUE)
  
  # Now, let's test a more complicated transformation:
  
  transformations <- "
  parameters: a,b,c,d,e,f,deltaA, deltaB, deltaC
  d = a + deltaA;
  e = b + deltaB;
  f = c + deltaC;
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
  
  
  ## Testing transformations against OpenMx
  modelLavaan <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + d*y6 + e*y7 + f*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
'
  
  fitLavaan <- sem(modelLavaan, 
                   data = PoliticalDemocracy, 
                   meanstructure = TRUE
  )
  
  transformations <- "
  parameters: a,b,c,d,e,f,deltaA, deltaB, deltaC
  d = a + deltaA;
  e = b + deltaB;
  f = c + deltaC;
  "
  rsemBFGS <- bfgs(lavaanModel = fitLavaan, 
                   modifyModel = modifyModel(transformations = transformations)
  )
  
  
  library(OpenMx)
  manifests <- c(paste0("x",1:3), paste0("y", 1:8))
  latents <- c("ind60", "dem60", "dem65")
  modelMx <- mxModel(
    mxPath(from = "ind60", 
           to = "dem60",
           labels = "dem60_ind60",
           free = c(TRUE)),
    mxPath(from = c("ind60",
                    "dem60"), 
           to = "dem65", 
           labels = c("dem65_ind60", "dem65_dem60"),
           free = TRUE),
    mxPath(from = "ind60", 
           to = paste0("x",1:3), 
           labels = paste0("l",1:3),
           free = c(FALSE, TRUE, TRUE), 
           values = c(1,1,1)),
    mxPath(from = "dem60", 
           to = paste0("y",1:4), 
           labels = c(NA, "a", "b", "c"),
           free = c(FALSE, TRUE, TRUE, TRUE), 
           values = c(1,1,1,1)),
    mxPath(from = "dem65", 
           to = paste0("y",5:8), 
           labels = c(NA, "d[1,1]", "e[1,1]", "f[1,1]"),
           free = c(FALSE, FALSE, FALSE, FALSE), 
           values = c(1,1,1,1)),
    mxAlgebra(expression = a+deltaA, name = "d"),
    mxAlgebra(expression = b+deltaB, name = "e"),
    mxAlgebra(expression = c+deltaC, name = "f"),
    mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 0, labels = "deltaA"),
    mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 0, labels = "deltaB"),
    mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 0, labels = "deltaC"),
    mxPath(from = latents, 
           to = latents, 
           labels = paste0(latents, "__", latents),
           arrows = 2, 
           free = TRUE),
    mxPath(from = manifests, 
           to = manifests, 
           labels = paste0(manifests, "__", manifests),
           arrows = 2,
           free = TRUE),
    mxData(observed = PoliticalDemocracy, type = "raw"),
    mxFitFunctionML(),
    latentVars = latents,
    manifestVars = manifests,
    mxPath(from = 'one', to = manifests),
    type = "RAM"
  )
  fitMx <- mxRun(modelMx)
  testthat::expect_equal(
    all(
      abs(
        omxGetParameters(fitMx)[c("a", "b", "c", "deltaA", "deltaB", "deltaC")] -
          unlist(rsemBFGS@parameters[c("a", "b", "c", "deltaA", "deltaB", "deltaC")])
      ) < .01), 
    TRUE)
  
  ## Let's test another transformation
  
  modelMx <- mxModel(
    mxPath(from = "ind60", 
           to = "dem60",
           labels = "dem60_ind60",
           free = c(TRUE)),
    mxPath(from = c("ind60",
                    "dem60"), 
           to = "dem65", 
           labels = c("dem65_ind60", "dem65_dem60"),
           free = TRUE),
    mxPath(from = "ind60", 
           to = paste0("x",1:3), 
           labels = paste0("l",1:3),
           free = c(FALSE, TRUE, TRUE), 
           values = c(1,1,1)),
    mxPath(from = "dem60", 
           to = paste0("y",1:4), 
           labels = c(NA, "a", "b", "c"),
           free = c(FALSE, TRUE, TRUE, TRUE), 
           values = c(1,1,1,1)),
    mxPath(from = "dem65", 
           to = paste0("y",5:8), 
           labels = c(NA, "d[1,1]", "e[1,1]", "f[1,1]"),
           free = c(FALSE, FALSE, FALSE, FALSE), 
           values = c(1,1,1,1)),
    mxAlgebra(expression = exp(dlog), name = "d"),
    mxAlgebra(expression = b+deltaB, name = "e"),
    mxAlgebra(expression = c+deltaC, name = "f"),
    mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 0, labels = "dlog"),
    mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 0, labels = "deltaB"),
    mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 0, labels = "deltaC"),
    mxPath(from = latents, 
           to = latents, 
           labels = paste0(latents, "__", latents),
           arrows = 2, 
           free = TRUE),
    mxPath(from = manifests, 
           to = manifests, 
           labels = paste0(manifests, "__", manifests),
           arrows = 2,
           free = TRUE),
    mxData(observed = PoliticalDemocracy, type = "raw"),
    mxFitFunctionML(),
    latentVars = latents,
    manifestVars = manifests,
    mxPath(from = 'one', to = manifests),
    type = "RAM"
  )
  fitMx <- mxRun(modelMx)
  omxGetParameters(fitMx)
  
  
  transformations <- "
  parameters: a,b,c,d,e,f,dlog, deltaB, deltaC
  d = exp(dlog);
  e = b + deltaB;
  f = c + deltaC;
  "
  rsemBFGS <- bfgs(lavaanModel = fitLavaan, 
                   modifyModel = modifyModel(transformations = transformations)
  )
  coef(rsemBFGS)
  
  testthat::expect_equal(
    all(
      abs(
        omxGetParameters(fitMx)[c("a", "b", "c", "dlog", "deltaB", "deltaC")] -
          unlist(rsemBFGS@parameters[c("a", "b", "c", "dlog", "deltaB", "deltaC")])
      ) < .01), 
    TRUE)
})
