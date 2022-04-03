test_that("scores work", {
  library(ipcr)
  library(lavaan)
  # test scores
  
  ## The following code is taken from the help function of iprc
  HS_data <- lavaan::HolzingerSwineford1939
  
  ## Remove observations with missing values
  HS_data <- HS_data[stats::complete.cases(HS_data), ]
  
  ## lavaan model syntac for a single group model
  m <- 'visual =~ x1 + x2 + x3
      textual =~ x4 + x5 + x6
      speed =~ x7 + x8 + x9'
  
  ## Fit the model
  model <- lavaan::cfa(model = m, data = HS_data)
  
  SEM <- SEMFromLavaan(model = model, rawData = HS_data)
  SEM <- fit(SEM)
  
  testthat::expect_equal(round(SEM$m2LL - (-2*as.numeric(logLik(model))),4),0)
  
  sc <- computeScores(SEM = SEM, individualFitfunction = individualLikelihoodRatioFit, raw = FALSE)
  sc2 <- get_scores(model)
  
  # Note: Scores are not independent of the fitting function; the fitting function used by 
  # ipcr is -.5 times the fitting function individualLikelihoodRatioFit
  
  testthat::expect_equal(any(abs(-.5*sc[,c("visual=~x2", "visual=~x3", "textual=~x5", "textual=~x6")] - 
                                   sc2[,c("visual..x2", "visual..x3", "textual..x5", "textual..x6")]) > .0001), FALSE)
  
  # compute ipcs
  hessian <- computeHessian(SEM = SEM, fitfunction = likelihoodRatioFit)
  # standard errors: sqrt(diag(solve(.5*hessian))); note: hessian is based in -2 log Likelihood; multiplying by .5 gives -log Likelihood
  bread_matrix <- ipcr:::bread_ipcr.lavaan(model)
  # standard errors: sqrt(diag(nrow(HS_data)^(-1)* bread_matrix))
  # -> nrow(HS_data)^(-1) bread_matrix = solve(.5*hessian)
  # -> bread_matrix = 2*nrow(HS_data)*hessian^(-1) = ((.5/N)*hessian)^(-1)
})
