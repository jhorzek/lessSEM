test_that("testing lsp", {
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  N <- 300
  l1 <- 1; l2 <- .9; l3 <- 0;
  v1 <- .2; v2 <- .3; v3 <- 1
  
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
  
  lavaanParameters <- lessSEM::getLavaanParameters(modelFit)
  lambdas = seq(0,.1,length.out = 3)
  thetas = seq(.1,.5,length.out = 2)
  
  rsemIsta <- lsp(lavaanModel = modelFit, 
                  regularized = paste0("f=~y",6:ncol(y)), 
                  lambdas = lambdas,
                  thetas = thetas,
                  control = controlIsta()
  )
  
  testthat::expect_equal(
    all(abs(rsemIsta@fits$m2LL[rsemIsta@fits$lambda == 0] - 
              (-2*logLik(modelFit))
    ) < 1e-2), 
    TRUE)
  
  # compare to smoothed version
  tuningParameters <- expand.grid(theta = thetas,
                                  lambda = lambdas
  )
  
  penaltyFunctionArguments <- list(
    eps = 1e-5,
    regularizedParameterLabels = paste0("f=~y",6:ncol(y))
  )
  
  #### Now we are ready to optimize! ####
  regsemApprox <- lessSEM:::.regularizeSEMWithCustomPenaltyRsolnp(lavaanModel = modelFit,
                                                                  individualPenaltyFunction = lessSEM:::.smoothLspValue,
                                                                  tuningParameters = tuningParameters,
                                                                  penaltyFunctionArguments = penaltyFunctionArguments)
  for(i in length(rsemIsta@fits$theta)){
    th <- rsemIsta@fits$theta[i]
    la <- rsemIsta@fits$lambda[i]
    
    testthat::expect_equal(
      abs(rsemIsta@fits$regM2LL[rsemIsta@fits$theta == th & 
                                  rsemIsta@fits$lambda == la] -
            regsemApprox@fits$regM2LL[regsemApprox@fits$theta == th & 
                                        regsemApprox@fits$lambda == la]
      )/rsemIsta@fits$regM2LL[rsemIsta@fits$theta == th & 
                                rsemIsta@fits$lambda == la] < .01,TRUE)
    testthat::expect_equal(
      all(
        abs(
          rsemIsta@parameters[rsemIsta@fits$theta == th &
                                rsemIsta@fits$lambda == la,rsemIsta@parameterLabels]-
            regsemApprox@parameters[regsemApprox@fits$theta == th &
                                      regsemApprox@fits$lambda == la,rsemIsta@parameterLabels]
        ) < .1
      ),
      TRUE
    )
  }
})
