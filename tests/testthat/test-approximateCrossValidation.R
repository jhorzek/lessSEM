test_that("testing approximate cross validation", {
  
  library(lavaan)
  library(regsem)
  library(lessSEM)
  #### Test lavaan ####
  # put variables on same scale for regsem
  HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
  mod <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
  outt = cfa(mod, HS, meanstructure = TRUE)
  
  aLOOCV <- aCV4lavaan(lavaanModel = outt, k = nrow(HS), raw = FALSE)
  
  exactLOOCV <- rep(NA, nrow(HS))
  for(i in 1:nrow(HS)){
    outt = cfa(mod, HS[-i,], meanstructure = TRUE)
    exactLOOCV[i] <- lessSEM:::computeIndividualM2LL(nObservedVariables = ncol(HS), 
                                                     rawData = as.numeric(HS[i,]),
                                                     impliedMeans = outt@implied$mean[[1]], 
                                                     impliedCovariance = outt@implied$cov[[1]])
  }
  
  testthat::expect_equal(max(abs(aLOOCV$leaveOutFits -exactLOOCV)) < .1, TRUE)
  
  plot(exactLOOCV, exactLOOCV, type = "l",
       xlab = "exact loocv", ylab = "approximated loocv")
  points(exactLOOCV, aLOOCV$leaveOutFits, col = "red")
  
})
