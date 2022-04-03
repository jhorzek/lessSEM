fit <- function(SEM){
  
  return(SEM$fit())
  
}


minus2LogLikelihood <- function(par, SEM, raw){
  
  if(any(getParameters(SEM)[names(par)] != par)) setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  SEM$fit()
  if(any(class(SEM) == "try-error") || !is.numeric(SEM$m2LL)){ return( 999999999 )}
  
  return(SEM$m2LL)
}

individualMinus2LogLikelihood <- function(par, SEM, data, raw){
  if(any(names(data) != SEM$manifestNames)) stop("SEM$manifestNames and colnames of data do not match!")
  if(any(getParameters(SEM, raw = raw)[names(par)] != par)) setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
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
  if(any(getParameters(SEM)[names(par)] != par)) setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  SEM$implied()
  
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

likelihoodRatioFitRegsem <- function(par, SEM, raw){
  N <- nrow(SEM$data$rawData)
  return(.5*N^(-1)*likelihoodRatioFit(par = par, SEM = SEM, raw = raw))
}

individualLikelihoodRatioFit <- function(par, SEM, data, raw){
  if(any(names(data) != SEM$manifestNames)) stop("Colnames of SEM$data$rawData and data do not match!")
  if(anyNA(data)) stop("Requires observed means and covariances for SEM object")
  m2LL <- individualMinus2LogLikelihood(par = par, SEM = SEM, data = data, raw = raw)
  isMissing <- is.na(data)
  saturatedFit <- computeIndividualM2LL(nObservedVariables = sum(!isMissing), 
                                        rawData = data[!isMissing],
                                        expectedMeans = SEM$data$obsMeans[!isMissing],
                                        expectedCovariance = SEM$data$obsCov[!isMissing,!isMissing])
  return(m2LL - saturatedFit)
}

minus2LogLikelihoodLASSO <- function(par, SEM, regularizedParameterLabels, lambda = NULL, lambda_ = NULL, eps){
  if(is.null(lambda) && is.null(lambda_)) stop("provide lambda or lambda_")
  if(is.null(lambda)) lambda <- lambda_
  m2LL <- minus2LogLikelihood(par, SEM)
  penalty <- lassoPenalty(par, regularizedParameterLabels, lambda, eps)
  return(m2LL + penalty)
}

likelihoodRatioFitLASSO <- function(par, SEM, regularizedParameterLabels, lambda, raw, eps){
  setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  tryFit <- try(fit(SEM), silent = TRUE)
  if(any(class(tryFit) == "try-error") || !is.numeric(SEM$m2LL)){ return( 999999999 )}
  
  penalty <- lassoPenalty(par, regularizedParameterLabels, lambda, eps)
  
  return(likelihoodRatioFit(par = par, SEM = SEM) + penalty)
}

likelihoodRatioFitRegsemLASSO <- function(par, SEM, regularizedParameterLabels, lambda = NULL, lambda_ = NULL, eps){
  if(is.null(lambda) && is.null(lambda_)) stop("provide lambda or lambda_")
  if(is.null(lambda)) lambda <- lambda_
  penalty <- lassoPenalty(par, regularizedParameterLabels, lambda, eps)
  
  return(likelihoodRatioFitRegsem(par = par, SEM = SEM) + .5*penalty)
}

individualLikelihoodRatioFitRegsemLASSO <- function(par, SEM, data, regularizedParameterLabels, lambda = NULL, lambda_ = NULL, raw, eps){
  if(is.null(lambda) && is.null(lambda_)) stop("provide lambda or lambda_")
  if(is.null(lambda)) lambda <- lambda_
  N <- nrow(SEM$data$rawData)
  penalty <- lassoPenalty(par, regularizedParameterLabels, lambda, eps)
  
  return(N^(-1)*(.5*individualLikelihoodRatioFit(par = par, SEM = SEM, data = data, raw = raw) + .5*penalty))
}

individualLikelihoodRatioFitRegsem <- function(par, SEM, data){
  
  N <- nrow(SEM$data$rawData)
  
  return(N^(-1)*(.5*individualLikelihoodRatioFit(par = par, SEM = SEM, data = data)))
}

lassoPenalty <- function(par, regularizedParameterLabels, lambda, eps){
  return(
    lambda*sum(sqrt(as.numeric(par)[names(par) %in% regularizedParameterLabels]^2 + eps))
  )
}
