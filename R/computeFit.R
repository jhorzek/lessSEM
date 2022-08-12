#' .fit
#' 
#' fits an object of class Rcpp_SEMCpp.
#' 
#' @param SEM model of class Rcpp_SEMCpp. 
#' @returns fitted SEM
.fit <- function(SEM){
  SEM$fit()
  return(SEM)
  
}

#' .fitFunction
#' 
#' internal function which returns the -2 log Likelihood of an object of class Rcpp_SEMCpp. This function can be used in optimizers
#' 
#' @param par labeled vector with parameter values
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of lessSEM is used.
#' @returns -2log-Likelihood
.fitFunction <- function(par, SEM, raw){
  SEM <- .setParameters(SEM = SEM, names(par), values = par, raw = raw)
  tryFit <- try(SEM$fit(), silent = TRUE)
  if(is(tryFit, "try-error") || is.na(SEM$m2LL) || any(eigen(SEM$S, only.values = TRUE)$values < 0)) {return(9999999999999)}
  return(SEM$m2LL)
}

#' .gradientFunction
#' 
#' internal function which returns the gradients of an object of class Rcpp_SEMCpp. This function can be used in optimizers
#' 
#' @param par labeled vector with parameter values
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of lessSEM is used.
#' @returns gradients of the model
.gradientFunction <- function(par, SEM, raw){
  failureReturns <- rep(999999999, length(par))
  names(failureReturns) <- names(par)
  
  SEM <- .setParameters(SEM = SEM, names(par), values = par, raw = raw)
  tryFit <- try(SEM$fit(), silent = TRUE)
  if(any(class(tryFit) == "try-error") || is.na(SEM$m2LL)) {return(failureReturns)}
  tryGradients <- try(.getGradients(SEM = SEM, raw = raw), silent = TRUE)
  if(any(class(tryGradients) == "try-error") || anyNA(tryGradients)) {return(failureReturns)}
  return(tryGradients)
}

#' .individualMinus2LogLikelihood
#' 
#' internal function which returns the -2 log Likelihood for a single subject
#' 
#' @param par labeled vector with parameter values
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param data vector with data points for this single individual
#' @param raw controls if the internal transformations of lessSEM is used.
#' @returns -2 log Likelihood for each subject in the data set
.individualMinus2LogLikelihood <- function(par, SEM, data, raw){
  if(any(names(data) != SEM$manifestNames)) stop("SEM$manifestNames and colnames of data do not match!")
  if(any(.getParameters(SEM, raw = raw)[names(par)] != par)) SEM <- .setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  if(anyNA(SEM$impliedCovariance) || anyNA(SEM$impliedMeans)) SEM$implied()
  isMissing <- is.na(data)
  return(
    computeIndividualM2LL(nObservedVariables = sum(!isMissing),  
                          rawData = data[!isMissing], 
                          impliedMeans = SEM$impliedMeans[!isMissing], 
                          impliedCovariance = SEM$impliedCovariance[!isMissing,!isMissing])
  )
  
}

#' .likelihoodRatioFit
#' 
#' internal function which returns the likelihood ratio fit statistic
#' 
#' @param par labeled vector with parameter values
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of lessSEM is used.
#' @returns likelihood ratio fit statistic
.likelihoodRatioFit <- function(par, SEM, raw){
  if(anyNA(SEM$rawData)) stop("likelihoodRatioFit currently only implemented for data without missings")
  if(any(.getParameters(SEM)[names(par)] != par)) SEM <- .setParameters(SEM, labels = names(par), values = as.numeric(par), raw = raw)
  SEM$fit()
  
  # compute fit
  N <- nrow(SEM$rawData)
  obsCov <- ((N-1)/N)*stats::cov(SEM$rawData)
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

