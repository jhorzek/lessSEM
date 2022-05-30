#' quasiNewtonBFGS
#' 
#' Performs quasi-Newton optimization with aCV4SEM:::BFGS approximated Hessian
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param individualPenaltyFunction penalty function which takes the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' returns a single value - the value of the penalty function for a single person. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionGradient gradients of the penalty function. Function should take the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' return a vector of the same length as parameters. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param individualPenaltyFunctionHessian Hessian of the penalty function. Function should take the current parameter values as first argument, the tuning parameters as second, and the penaltyFunctionArguments as third argument and 
#' return a matrix with (length as parameters)^2 number of elements. If the true penalty function is non-differentiable (e.g., lasso) a smooth
#' approximation of this function should be provided.
#' @param currentTuningParameters data.frame with tuning parameter values. 
#' @param penaltyFunctionArguments arguments passed to individualPenaltyFunction, individualPenaltyFunctionGradient, and individualPenaltyFunctionHessian
#' @param initialHessian option to provide an initial Hessian to the optimizer
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param epsOut Stopping criterion for outer iterations
#' @param epsIn Stopping criterion for inner iterations
#' @param convergenceCriterion which convergence criterion should be used for the outer iterations? possible are "GLMNET", "gradients", "fitChange"
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
quasiNewtonBFGS <- function(SEM, 
                            individualPenaltyFunction, 
                            individualPenaltyFunctionGradient,
                            individualPenaltyFunctionHessian,
                            currentTuningParameters,
                            penaltyFunctionArguments,
                            initialHessian,
                            stepSize,
                            sig,
                            gam,
                            maxIterOut,
                            maxIterIn,
                            maxIterLine,
                            epsOut,
                            epsIn,
                            convergenceCriterion,
                            verbose = 0){
  
  # Setup
  N <- nrow(SEM$rawData)
  
  # save current state
  initialParameters <- aCV4SEM:::getParameters(SEM, raw = TRUE)
  if(is.null(initialHessian)){
    initialHessian <- aCV4SEM:::getHessian(SEM = SEM, raw = TRUE)
  }
  initialGradients <- try(aCV4SEM:::getGradients(SEM = SEM, raw = TRUE), silent = TRUE)
  if(is(initialGradients, "try-error")){stop("Could not compute gradients at initial parameter values.")}
  
  ## check penalty functions
  penaltyValue <- N*individualPenaltyFunction(initialParameters, currentTuningParameters, penaltyFunctionArguments)
  if(length(penaltyValue) != 1 || !is.numeric(penaltyValue)) stop("individualPenaltyFunctionGradient must return a numeric scalar.")
  penaltyGradient <- N*individualPenaltyFunctionGradient(initialParameters, currentTuningParameters, penaltyFunctionArguments)
  if(length(penaltyGradient) != length(initialParameters) || any(!names(penaltyGradient) %in% names(initialGradients))) stop("individualPenaltyFunctionGradient must return a labeled vector with the same names as the parameters")
  penaltyHessian <- N*individualPenaltyFunctionHessian(initialParameters, currentTuningParameters, penaltyFunctionArguments)
  if(nrow(penaltyHessian) != length(initialParameters) || 
     ncol(penaltyHessian) != length(initialParameters) || 
     any(!rownames(penaltyHessian) %in% names(initialGradients)) || 
     any(!colnames(penaltyHessian) %in% names(initialGradients))) stop("individualPenaltyFunctionHessian must return a symmetric matrix with the parameter names as rownames and colnames.")
  
  # add penalty
  initialGradients <- initialGradients + penaltyGradient
  initialHessian <- initialHessian + penaltyHessian
  
  # initialize parameter values
  newParameters <- initialParameters
  parameterLabels <- names(newParameters)
  
  #### start optimizing ####
  
  # outer loop: optimize parameters
  iterOut <- 0
  
  converged <- TRUE
  
  # initialize parameters
  newParameters <- initialParameters
  # initialize initial gradients
  newGradients <- initialGradients
  # initialize Hessian
  newHessian <- initialHessian
  eigenDecomp <- eigen(newHessian)
  if(any(eigenDecomp$values < 0)){
    warning("Initial Hessian is not positive definite. Using identity matrix instead. If you are using the numDeriv approximation, this may be the cause of the non-positive-definiteness. This is especially common if the penalty function is still relatively non-smooth.")
    newHessian <- diag(1, nrow(initialHessian))
    rownames(newHessian) <- parameterLabels
    colnames(newHessian) <- parameterLabels
  }
  
  newM2LL <- SEM$m2LL
  newRegM2LL <- newM2LL + penaltyValue
  
  while(iterOut < maxIterOut){
    iterOut <- iterOut + 1
    
    if(iterOut > 1){
      parameterChange <- (newParameters - oldParameters)
    }
    
    oldParameters <- newParameters
    oldGradients <- newGradients
    oldHessian <- newHessian
    
    oldM2LL <- newM2LL
    if(iterOut > 1){
      regM2LLChange <- newRegM2LL - oldRegM2LL
    }
    oldRegM2LL <- newRegM2LL
    
    # outer breaking condition
    if(iterOut > 1){
      # requires a direction vector d; therefore, we cannot evaluate the
      # convergence in the first iteration
      # requires a direction vector d; therefore, we cannot evaluate the
      # convergence in the first iteration
      if(convergenceCriterion == "GLMNET") converged <- try(max(diag(diag(oldHessian))%*%direction^2) < epsOut,
                                                            silent = TRUE)
      if(convergenceCriterion == "gradients") converged <- try(all(abs(newGradients/N) < # gradient are based on likelihood and therefore 
                                                                     # increase with the number of subjects
                                                                     epsOut), 
                                                               silent = TRUE)
      if(convergenceCriterion == "fitChange") converged <- try(abs(regM2LLChange) < epsOut, 
                                                               silent = TRUE)
      
      if(is(converged, "try-error") || is.na(converged)){
        converged <- FALSE
        # save last working parameters
        newParameters <- oldParameters
        newGradients <- oldGradients
        newHessian <- oldHessian
        warning(paste("The model did NOT CONVERGE for lambda = ", lambda, "."))
        break
      }
      if(converged){
        break
      }
      
      if(verbose > 0){
        
        cat(paste0("\r",
                   "## [", iterOut,
                   "] m2LL: ", sprintf('%.4f',newM2LL),
                   " | regM2LL:  ", sprintf('%.4f',newRegM2LL),
                   " | regM2LL change:  ", sprintf('%.4f',regM2LLChange),
                   " | max. gradient:  ", sprintf('%.5f',max(abs(newGradients/N))),
                   " ##")
        )
        flush.console()
      }
    }
    
    
    ## inner loop: optimize directions
    direction <- c(-solve(newHessian)%*%newGradients)
    names(direction) <- parameterLabels
    
    # perform Line Search
    step <- quasiNewtonLineSearch(SEM = SEM, 
                                  N = N,
                             individualPenaltyFunction = individualPenaltyFunction,
                             currentTuningParameters = currentTuningParameters,
                             penaltyFunctionArguments = penaltyFunctionArguments,
                             oldParameters = oldParameters, 
                             oldM2LL = oldM2LL, 
                             oldGradients = oldGradients,
                             oldHessian = oldHessian,
                             direction = direction,
                             stepSize= stepSize, 
                             sig = sig,
                             gam = gam, 
                             maxIterLine = maxIterLine)
    
    # extract elements:
    newParameters <- step$newParameters
    newGradients <- step$newGradients + N*individualPenaltyFunctionGradient(newParameters, 
                                                                            currentTuningParameters,
                                                                          penaltyFunctionArguments)
    
    # update model: set parameter values and compute
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(newParameters), values = newParameters, raw = TRUE)
    SEM <- fit(SEM)
    
    # get fit
    newM2LL <- SEM$m2LL
    newRegM2LL <- newM2LL + N*individualPenaltyFunction(newParameters, 
                                                        currentTuningParameters,
                                                      penaltyFunctionArguments)
    
    # Approximate Hessian using aCV4SEM:::BFGS
    newHessian <- aCV4SEM:::BFGS(oldParameters = oldParameters, 
                                 oldGradients = oldGradients, 
                                 oldHessian = oldHessian, 
                                 newParameters = newParameters, 
                                 newGradients = newGradients)
  }
  # warnings
  if(iterOut == maxIterOut){
    warning(paste("For currentTuningParameters = ", paste0(currentTuningParameters, collapse = "; "), "the maximum number of iterations was reached."))
  }
  
  return(list("SEM" = SEM, 
              "parameters" = newParameters, 
              "m2LL" = newM2LL,
              "regM2LL" = newRegM2LL,
              "gradients" = newGradients, 
              "Hessian" = newHessian, 
              "convergence" = converged))
}




#### Line search ####

#' quasiNewtonLineSearch
#'
#' performs a line search procedure adapted from Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421 Equation 20.
#'
#'
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param N sample size
#' @param regularizedParameterLabels labels of regularized parameters
#' @param adaptiveLassoWeights vector with weights for adaptive LASSO. Set to NULL if not using adaptive LASSO
#' @param lambda lasso tuning parameter. Higher values = higher penalty
#' @param oldParameters parameters of previous iteration
#' @param oldM2LL -2 log likelihood of previous iteration
#' @param oldGradients gradients of previous iteration
#' @param oldHessian Hessian approximation of previous iteration
#' @param direction vector with step direction
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = theta_k + Stepsize \* Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterLine maximal number of iterations for line search
quasiNewtonLineSearch <- function(SEM, 
                                  N,
                                  individualPenaltyFunction,
                                  currentTuningParameters,
                                  penaltyFunctionArguments,
                                  oldParameters, 
                                  oldM2LL, 
                                  oldGradients,
                                  oldHessian,
                                  direction,
                                  stepSize, 
                                  sig,
                                  gam, 
                                  maxIterLine){
  newGradients <- NULL
  # get penalized M2LL for step size 0:
  pen_0 <- N*individualPenaltyFunction(oldParameters, 
                                       currentTuningParameters,
                                       penaltyFunctionArguments)
  f_0 <- oldM2LL + pen_0
  
  i <- 0
  stepSizeInit <- stepSize
  if(stepSizeInit >= 1) stepSizeInit <- .9
  # sometimes the optimizer gets stuck in one specific location;
  # here, it can help to try different step sizes
  stepSizeInit <- runif(1,.1,stepSizeInit)
  while(TRUE){
    stepSize <- stepSizeInit^i
    newParameters <- oldParameters+stepSize*direction
    
    # compute new fitfunction value
    newM2LL <- try(aCV4SEM:::fit(aCV4SEM:::setParameters(SEM, names(newParameters), newParameters, raw = TRUE))$m2LL, silent = TRUE)
    if(!is.finite(newM2LL)){
      i <- i+1
      next
    }
    if(is(newM2LL, "try-error")){
      stop("Error in line search: Could not compute fit at current location.")
    }
    
    
    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    p_new <- N*individualPenaltyFunction(newParameters, currentTuningParameters, penaltyFunctionArguments)
    f_new <- newM2LL + p_new
    
    # test line search criterion
    lineCriterion <- f_new - f_0 <= sig*stepSize*(t(oldGradients)%*%direction + gam*t(direction)%*%oldHessian%*%direction + p_new - pen_0)
    if(lineCriterion){
      # check if covariance matrix is positive definite
      if(any(eigen(SEM$S, only.values = TRUE)$values < 0)){
        i <- i+1
        next
      }
      if(any(eigen(SEM$impliedCovariance, only.values = TRUE)$values < 0)){
        i <- i+1
        next
      }
      # check if gradients can be computed at the new location; this can often cause issues
      newGradients <-  try(aCV4SEM:::getGradients(SEM, raw = TRUE), silent = TRUE)
      
      if(is(newGradients, "try-error")) {
        car("Was error\n")
        i <- i+1
        next
      }
      break
    }
    i <- i+1
    if(i >= maxIterLine){
      break
    }
  }
  if(is.null(newGradients) || is(newGradients, "try-error")){
    newGradients <- try(aCV4SEM:::getGradients(SEM, raw = TRUE), silent = TRUE)
  }
  return(
    list("stepSize" = stepSize,
         "newParameters" = newParameters,
         "newGradients" = newGradients)
  )
}
