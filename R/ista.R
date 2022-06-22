#' istaLasso
#' 
#' Optimize LASSO regularized SEM with ista optimizer 
#' @param SEM object of class SEMCpp
#' @param lambda tuning parameter lambda
#' @param weights weights given to penalty for each parameter
#' @param i_k currently not used
#' @param L0 controls the initial step size and change in step size (L0*(eta^k), where k is the iteration) 
#' @param eta controls change of step size in ista (L0*(eta^k), where k is the iteration) 
#' @param maxIterOut maximal number of outer iterations
#' @param maxIterIn maximal number of inner iterations
#' @param breakOuter change in fit below which convergence is assumend
#' @param raw should parameters be returned in raw format?
istaLasso <- function(SEM,
                      lambda,
                      weights,
                      i_k = 2,
                      L0 = 1,
                      eta = 2,
                      maxIterOut = 1000,
                      maxIterIn = 1000,
                      breakOuter = .0000001,
                      raw = FALSE){
  
  startingValues <- linr:::getParameters(SEM, raw = TRUE)
  if(any(! names(startingValues) %in% names(weights))) stop("weights vector must have the same names as the parameters!")
  weights <- weights[names(startingValues)]
  
  # optimize
  optPar <- SEM$LASSO(
    startingValues,
    lambda,
    weights,
    i_k,
    L0,
    eta,
    maxIterOut,
    maxIterIn,
    breakOuter
  )
  if(raw) return(optPar)
  
  SEM <-linr:::setParameters(SEM, labels = names(optPar), values = optPar, raw = TRUE)
  return(linr:::getParameters(SEM, raw = FALSE))
}