test_that("testing gradients", {
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  
  model1 <- ' 
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
    y3 ~~ d*y7
    y4 ~~ d*y8
    y6 ~~ y8
    
    ind60 ~ 1
    y1 ~ 0; y2 ~ 0; y3 ~ 0
    
'
  
  model <- sem(model1, data = PoliticalDemocracy, meanstructure = TRUE)
  scores <- -2*lavScores(model)
  gradients <- apply(scores,2,sum)
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model)
  SEM <- lessSEM:::.fit(SEM)
  
  gradientsAnalytic <- lessSEM:::.getGradients(SEM = SEM, raw = FALSE)

  testthat::expect_equal(round(sum(abs(gradientsAnalytic - gradients)),4),0)
  
  PoliticalDemocracyWithNA <- PoliticalDemocracy
  rows <- 1:nrow(PoliticalDemocracyWithNA)
  cols <- 1:ncol(PoliticalDemocracyWithNA)
  combinations <- expand.grid(rows,cols)
  
  set.seed(123)
  makeNA <- sample(1:nrow(combinations), 50)
  for(i in 1:length(makeNA)){
    PoliticalDemocracyWithNA[combinations[makeNA[i],1], combinations[makeNA[i],2]] <- NA
  }
  
  model <- sem(model1, data = PoliticalDemocracyWithNA, meanstructure = TRUE, missing = "ML")
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model)
  SEM <- lessSEM:::.fit(SEM)
  
  scores <- -2*lavScores(model)
  gradients <- apply(scores,2,sum)
  
  gradientsAnalytic <- lessSEM:::.getGradients(SEM, raw = FALSE)
  testthat::expect_equal(round(sum(abs(gradientsAnalytic - gradients)),4),0)
})
