test_that("testing scores", {
  testthat::skip_on_cran()
  library(lavaan)
  library(lessSEM)
  set.seed(123)

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
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
  
  model <- sem(model1, data = PoliticalDemocracy, meanstructure = TRUE)
  lavaanScores <- -2*lavaan::lavScores(model)
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model)
  SEM <- lessSEM:::.fit(SEM)
  
  scoresAnalytic <- lessSEM:::.getScores(SEM, raw = FALSE)

  testthat::expect_equal(round(sum(abs(scoresAnalytic - lavaanScores)),4),0)
  
  PoliticalDemocracyWithNA <- PoliticalDemocracy
  
  rows <- 1:nrow(PoliticalDemocracyWithNA)
  cols <- 1:ncol(PoliticalDemocracyWithNA)
  combinations <- expand.grid(rows,cols)
  
  makeNA <- sample(1:nrow(combinations), 50)
  for(i in 1:length(makeNA)){
    PoliticalDemocracyWithNA[combinations[makeNA[i],1], combinations[makeNA[i],2]] <- NA
  }
  
  model <- sem(model1, data = PoliticalDemocracyWithNA, meanstructure = TRUE, missing = "ML")
  lavaanScores <- -2*lavaan::lavScores(model)
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = model)
  SEM <- lessSEM:::.fit(SEM)
  
  scoresAnalytic <- lessSEM:::.getScores(SEM, raw = FALSE)
  testthat::expect_equal(round(sum(abs(scoresAnalytic - lavaanScores)),4),0)
})
