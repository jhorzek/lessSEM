#' .SEMdata
#' 
#' internal function. Creates internal data representation
#' @param rawData matrix with raw data set
#' @returns list with internal representation of data
#' @keywords internal
.SEMdata <- function(rawData){
  
  # extract unique missingness patterns
  
  isMissing <- is.na(rawData)
  uniqueMissingPatterns <- unique(isMissing)
  individualMissingPatternID <- rep(NA, nrow(rawData))
  
  missingSubsets <- vector("list", nrow(uniqueMissingPatterns))
  
  for(mrow in 1:nrow(uniqueMissingPatterns)){
    
    pattern <- uniqueMissingPatterns[mrow, ]
    
    individuals <- logicalMatch(isMissing, pattern)
    
    individualMissingPatternID[individuals] <- mrow
    
    # build subsets
    if(length(individuals) > 1){
      missingSubsets[[mrow]]$N <- length(individuals)
      missingSubsets[[mrow]]$observed <- sum(!pattern)
      missingSubsets[[mrow]]$notMissing <- which(!pattern)-1
      missingSubsets[[mrow]]$covariance <- ((length(individuals)-1)/length(individuals))*stats::cov(rawData[individuals,!pattern,drop = FALSE])
      missingSubsets[[mrow]]$means <- apply(rawData[individuals,!pattern,drop = FALSE], 2, mean)
      missingSubsets[[mrow]]$rawData <- rawData[individuals,,drop=FALSE]
      missingSubsets[[mrow]]$objectiveValue <- NA
    }else{
      missingSubsets[[mrow]]$N <- length(individuals)
      missingSubsets[[mrow]]$observed <- sum(!pattern)
      missingSubsets[[mrow]]$notMissing <- which(!pattern)-1
      missingSubsets[[mrow]]$covariance <- matrix(NA, nrow = sum(!pattern), ncol = sum(!pattern))
      missingSubsets[[mrow]]$means <-  matrix(NA, nrow = sum(!pattern), ncol = 1)
      missingSubsets[[mrow]]$rawData <- rawData[individuals,,drop=FALSE]
      missingSubsets[[mrow]]$objectiveValue <- NA
    }
  }
  
  dataList <- list("uniqueMissingPatterns" = uniqueMissingPatterns,
                   "individualMissingPatternID" = individualMissingPatternID,
                   "missingSubsets" = missingSubsets
  )
  
  return(dataList)
  
}

#' .SEMdataWLS
#' 
#' internal function. Creates internal data representation
#' @param rawData matrix with raw data set
#' @param lavaanModel lavaan model
#' @returns list with internal representation of data
#' @keywords internal
.SEMdataWLS <- function(rawData, lavaanModel){
  
  # unique missingness patterns -> none in case of WLS, but we will 
  # use the same data structure
  N <- lavInspect(lavaanModel, "nobs")
  observedMean <- lavInspect(lavaanModel, "sampstat")$mean
  observedCov <- lavInspect(lavaanModel, "sampstat")$cov
  if(is.null(observedMean)) 
    observedMean <- matrix(NA, nrow = 1, ncol = ncol(observedCov))
  
  uniqueMissingPatterns <- matrix(FALSE, nrow = 1, ncol = ncol(observedCov))
  individualMissingPatternID <- rep(1, N)
  
  missingSubsets <- vector("list", nrow(uniqueMissingPatterns))
  
  # build subsets
  mrow <- 1
  
  missingSubsets[[mrow]]$N <- N
  missingSubsets[[mrow]]$observed <- sum(!uniqueMissingPatterns[mrow,])
  missingSubsets[[mrow]]$notMissing <- which(!uniqueMissingPatterns[mrow,])-1
  missingSubsets[[mrow]]$covariance <- observedCov
  missingSubsets[[mrow]]$means <- observedMean
  missingSubsets[[mrow]]$rawData <- rawData
  missingSubsets[[mrow]]$objectiveValue <- NA
  
  dataList <- list("uniqueMissingPatterns" = uniqueMissingPatterns,
                   "individualMissingPatternID" = individualMissingPatternID,
                   "missingSubsets" = missingSubsets
  )
  
  return(dataList)
  
}

