test_that("testing mcp", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("lslx")
  library(lslx)
  library(lavaan)
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
  
  # fit model using lslx
  lslxModelSyntax <- paste0(paste0('fix(1)*', yNames[1], ' + ', paste0(yNames[2:5], collapse = " + "), " <=: f"),"\n",
                            paste0(paste0(yNames[6:length(yNames)], collapse = " + "), " <~: f"),"\n",
                            paste0(yNames, collapse = " + "), " <= 1"
  )
  fitLslx <- lslx$new(model = lslxModelSyntax,
                      sample_cov = stats::cov(y),
                      sample_size = nrow(y)
  )
  
  #fitLslx$penalize_coefficient(name = paste0("y", 6:ncol(y)," <- f"))
  
  lambdas <- seq(0,.3,length.out = 5)
  thetas <- seq(2.1,3,length.out = 2)
  fitLslx$fit(penalty_method = "mcp",
              lambda_grid = lambdas, 
              delta_grid = thetas,
              loss = "ml")
  
  # extract fits
  lslxParameter <- expand.grid(
    theta = thetas,
    lambda = lambdas
  )
  lslxParameter <- cbind(
    lslxParameter,
    matrix(NA, 
           nrow = nrow(lslxParameter),
           ncol = length(fitLslx$extract_coefficient(lambda = 0, 
                                                     delta = thetas[1])),
           dimnames = list(NULL, names(fitLslx$extract_coefficient(lambda = 0, 
                                                                   delta = thetas[1]))))
  )
  
  fits <- expand.grid(
    theta = thetas,
    lambda = lambdas,
    m2LL = NA,
    AIC = NA,
    BIC = NA,
    penalty = NA
  )
  
  regularized <- paste0("y", 6:ncol(y),"<-f/g") 
  
  for(l in 1:length(lambdas)){
    for(th in 1:length(thetas)){
      pars <- fitLslx$extract_coefficient(lambda = lambdas[l],
                                          delta = thetas[th])
      
      lslxParameter[lslxParameter$lambda == lambdas[l] &
                      lslxParameter$theta == thetas[th]
                    ,names(pars)] <- pars
      
      # lslx seems to not compute the implied means correctly here:
      #impl_means <- fitLslx$extract_implied_mean(lambda = lambdas[l])$g
      impl_means <- apply(y,2,mean)
      impl_cov <- fitLslx$extract_implied_cov(lambda = lambdas[l],
                                              delta = thetas[th])$g
      coefficents <- fitLslx$extract_coefficient(lambda = lambdas[l],
                                                 delta = thetas[th])
      regularizedCoefficients <- coefficents[regularized]
      
      m2LL <- -2*sum(mvtnorm::dmvnorm(x = y, mean = impl_means, sigma = impl_cov, log = TRUE))
      
      fits$m2LL[
        lslxParameter$lambda == lambdas[l] &
          lslxParameter$theta == thetas[th]
      ] <- m2LL
      
      fits$AIC[
        lslxParameter$lambda == lambdas[l] &
          lslxParameter$theta == thetas[th]
      ] <- m2LL + 2*(length(coefficents) - sum(regularizedCoefficients == 0) -1) # -1 because lslx also reports the fixed loading
      
      fits$BIC[
        lslxParameter$lambda == lambdas[l] &
          lslxParameter$theta == thetas[th]
      ] <-  m2LL + log(N)*(length(coefficents) - sum(regularizedCoefficients == 0) -1)
      
      
      fits$penalty[
        lslxParameter$lambda == lambdas[l] &
          lslxParameter$theta == thetas[th]
      ] <- 0
      for(r in regularized){
        fits$penalty[
          lslxParameter$lambda == lambdas[l] &
            lslxParameter$theta == thetas[th]
        ] <- fits$penalty[
          lslxParameter$lambda == lambdas[l] &
            lslxParameter$theta == thetas[th]
        ] + N * mcpPenalty_C(coefficents[r], lambda_p = lambdas[l], theta = thetas[th])
      }
    }
  }
  
  label <- lavaan2lslxLabels(lavaanModel = modelFit)
  label$lavaanLabels <- label$lavaanLabels[!grepl(pattern = "<-1/g", label$lslxLabels)]
  label$lslxLabels <- label$lslxLabels[!grepl(pattern = "<-1/g", label$lslxLabels)]
  
  # replicate with regularizedSEM
  lavaanParameters <- lessSEM::getLavaanParameters(modelFit)
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  
  rsemIsta <- mcp(lavaanModel = modelFit, 
                  regularized = regularizedLavaan, 
                  thetas = thetas,
                  lambdas = lambdas,
                  control = controlIsta(convCritInner = 0, verbose = 0)
  )
  
  testthat::expect_equal(
    all(round(rsemIsta@fits$m2LL[rsemIsta@fits$lambda == 0] -
                -2*logLik(modelFit)
              ,4)==0),TRUE)
  
  for(i in 1:length(rsemIsta@fits$theta)){
    th <- rsemIsta@fits$theta[i]
    la <- rsemIsta@fits$lambda[i]
    testthat::expect_equal(
      round(rsemIsta@fits$regM2LL[rsemIsta@fits$theta == th &
                                    rsemIsta@fits$lambda == la] -
              (fits$m2LL[fits$theta == th &
                           fits$lambda == la]
               + fits$penalty[fits$theta == th &
                                fits$lambda == la]
              )
      )==0,TRUE)
    
    testthat::expect_equal(
      max(abs(lslxParameter[lslxParameter[,"theta"] == th & 
                    lslxParameter[,"lambda"] == la,label$lslxLabels] -
      rsemIsta@parameters[
        rsemIsta@parameters$lambda == la &  rsemIsta@parameters$theta == th
        ,label$lavaanLabels])) < 1e-1, TRUE)
    
  }
  
})
