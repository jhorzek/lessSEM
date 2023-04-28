#' .fit
#' 
#' fits an object of class Rcpp_SEMCpp.
#' 
#' @param SEM model of class Rcpp_SEMCpp. 
#' @returns fitted SEM
#' @keywords internal
.fit <- function(SEM){
  SEM$fit()
  return(SEM)
}

#' .fitFunction
#' 
#' internal function which returns the objective value of the fitting function of an object of class Rcpp_SEMCpp. This function can be used in optimizers
#' 
#' @param par labeled vector with parameter values
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of lessSEM is used.
#' @returns objective value of the fitting function
#' @keywords internal
.fitFunction <- function(par, SEM, raw){
  SEM <- .setParameters(SEM = SEM, names(par), values = par, raw = raw)
  tryFit <- try(SEM$fit(), silent = TRUE)
  if(is(tryFit, "try-error") || is.na(SEM$objectiveValue) || any(eigen(SEM$S, only.values = TRUE)$values < 0)) {return(9999999999999)}
  return(SEM$objectiveValue)
}

#' .gradientFunction
#' 
#' internal function which returns the gradients of an object of class Rcpp_SEMCpp. This function can be used in optimizers
#' 
#' @param par labeled vector with parameter values
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of lessSEM is used.
#' @returns gradients of the model
#' @keywords internal
.gradientFunction <- function(par, SEM, raw){
  failureReturns <- rep(999999999, length(par))
  names(failureReturns) <- names(par)
  
  SEM <- .setParameters(SEM = SEM, names(par), values = par, raw = raw)
  tryFit <- try(SEM$fit(), silent = TRUE)
  if(any(class(tryFit) == "try-error") || is.na(SEM$objectiveValue)) {return(failureReturns)}
  tryGradients <- try(.getGradients(SEM = SEM, raw = raw), silent = TRUE)
  if(any(class(tryGradients) == "try-error") || anyNA(tryGradients)) {return(failureReturns)}
  return(tryGradients)
}


#' .standardErrors
#' 
#' compute the standard errors of a fitted SEM. IMPORTANT: Assumes that the 
#' SEM has been fitted and the parameter estimates are at the ordinary maximum
#' likelihood estimates
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of lessSEM is used. If 
#' set to TRUE, the standard errors will be returned for the internally used 
#' parameter specification
#' @return a vector with standard errors
#' @keywords internal
.standardErrors <- function(SEM, raw){
  
  Hessian <- .getHessian(SEM = SEM, raw = raw)
  
  # Note: We are minimizing the 2 times negative log likelihood.
  # The Fisher Information is the negative expected Hessian of the log likelihood
  # The Hessian of the 2 times negative log likelihood is therefore 2*"observed Fisher Information"
  # and 2 times it's inverse is the covariance matrix of the parameters
  nFisherInformation <- .5*(Hessian)
  standardErrors <- sqrt(diag(solve(nFisherInformation)))
  
  return(standardErrors)
}

#' .likelihoodRatioFit
#' 
#' internal function which returns the likelihood ratio fit statistic
#' 
#' @param par labeled vector with parameter values
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param raw controls if the internal transformations of lessSEM is used.
#' @returns likelihood ratio fit statistic
#' @keywords internal
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

