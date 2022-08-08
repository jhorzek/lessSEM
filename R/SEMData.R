#' .SEMdata
#' 
#' internal function. Creates internal data representation
#' @param rawData matrix with raw data set
#' @returns list with internal representation of data
.SEMdata <- function(rawData){
  
  # extract unique missingness patterns
  
  isMissing <- is.na(rawData)
  uniqueMissingPatterns <- unique(isMissing)
  individualMissingPatternID <- rep(NA, nrow(rawData))
  
  missingSubsets <- vector("list", nrow(uniqueMissingPatterns))
  
  for(mrow in 1:nrow(uniqueMissingPatterns)){
    individuals <- which(apply(isMissing, 1, function(x) sum(abs(x-uniqueMissingPatterns[mrow,]))) == 0)
    individualMissingPatternID[individuals] <- mrow
    
    # build subsets
    if(length(individuals) > 1){
      missingSubsets[[mrow]]$N <- length(individuals)
      missingSubsets[[mrow]]$observed <- sum(!uniqueMissingPatterns[mrow,])
      missingSubsets[[mrow]]$notMissing <- which(!uniqueMissingPatterns[mrow,])-1
      missingSubsets[[mrow]]$covariance <- ((length(individuals)-1)/length(individuals))*stats::cov(rawData[individuals,!uniqueMissingPatterns[mrow,]])
      missingSubsets[[mrow]]$means <- apply(rawData[individuals,!uniqueMissingPatterns[mrow,]], 2, mean)
      missingSubsets[[mrow]]$rawData <- rawData[individuals,,drop=FALSE]
      missingSubsets[[mrow]]$m2LL <- NA
    }else{
      missingSubsets[[mrow]]$N <- length(individuals)
      missingSubsets[[mrow]]$observed <- sum(!uniqueMissingPatterns[mrow,])
      missingSubsets[[mrow]]$notMissing <- which(!uniqueMissingPatterns[mrow,])-1
      missingSubsets[[mrow]]$covariance <- matrix(NA, nrow = sum(!uniqueMissingPatterns[mrow,]), ncol = sum(!uniqueMissingPatterns[mrow,]))
      missingSubsets[[mrow]]$means <-  matrix(NA, nrow = sum(!uniqueMissingPatterns[mrow,]), ncol = 1)
      missingSubsets[[mrow]]$rawData <- rawData[individuals,,drop=FALSE]
      missingSubsets[[mrow]]$m2LL <- NA
    }
  }
  
  dataList <- list("uniqueMissingPatterns" = uniqueMissingPatterns,
                   "individualMissingPatternID" = individualMissingPatternID,
                   "missingSubsets" = missingSubsets
  )
  
  return(dataList)
  
}

