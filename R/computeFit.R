fit <- function(SEM){
  SEM$fit()
  return(SEM)
  
}

fitFunction <- function(par, SEM, raw){
  SEM <- setParameters(SEM = SEM, names(par), values = par, raw = raw)
  tryFit <- try(SEM$fit(), silent = TRUE)
  if(any(class(tryFit) == "try-error") || is.na(SEM$m2LL)) {return(9999999999999)}
  return(SEM$m2LL)
}

derivativeFunction <- function(par, SEM, raw){
  failureReturns <- rep(999999999, length(par))
  names(failureReturns) <- names(par)
  
  SEM <- setParameters(SEM = SEM, names(par), values = par, raw = raw)
  tryFit <- try(SEM$fit(), silent = TRUE)
  if(any(class(tryFit) == "try-error") || is.na(SEM$m2LL)) {return(failureReturns)}
  tryGradients <- try(getGradients(SEM = SEM, raw = raw), silent = TRUE)
  if(any(class(tryGradients) == "try-error") || anyNA(tryGradients)) {return(failureReturns)}
  return(tryGradients)
}

minus2LogLikelihood <- function(par, SEM, raw){
  
  if(any(getParameters(SEM)[names(par)] != par)) setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  SEM$fit()
  if(any(class(SEM) == "try-error") || !is.numeric(SEM$m2LL)){ return( 999999999 )}
  
  return(SEM$m2LL)
}

individualMinus2LogLikelihood <- function(par, SEM, data, raw){
  if(any(names(data) != SEM$manifestNames)) stop("SEM$manifestNames and colnames of data do not match!")
  if(any(getParameters(SEM, raw = raw)[names(par)] != par)) SEM <- setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  if(anyNA(SEM$impliedCovariance) || anyNA(SEM$impliedMeans)) SEM$implied()
  isMissing <- is.na(data)
  return(
    computeIndividualM2LL(nObservedVariables = sum(!isMissing),  
                          rawData = data[!isMissing], 
                          impliedMeans = SEM$impliedMeans[!isMissing], 
                          impliedCovariance = SEM$impliedCovariance[!isMissing,!isMissing])
  )
  
}

likelihoodRatioFit <- function(par, SEM, raw){
  if(anyNA(SEM$rawData)) stop("likelihoodRatioFit currently only implemented for data without missings")
  if(any(getParameters(SEM)[names(par)] != par)) SEM <- setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  SEM$fit()
  
  # compute fit
  N <- nrow(SEM$rawData)
  obsCov <- ((N-1)/N)*cov(SEM$rawData)
  obsMeans <- apply(SEM$rawData,2,mean)
  expectedCovInverse <- solve(SEM$impliedCovariance)
  lrValue <- N*(log(det(SEM$impliedCovariance)) + 
                  sum(diag(obsCov%*%expectedCovInverse)) + 
                  matrix(SEM$impliedMeans - obsMeans, nrow = 1)%*%expectedCovInverse%*%matrix(SEM$impliedMeans - obsMeans, ncol = 1) -
                  log(det(obsCov)) -
                  length(obsMeans)
  )
  return(lrValue)
}

