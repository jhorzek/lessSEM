#' GLMNET
#' 
#' Performs GLMNET (see Friedman, 2010 & Yuan, 2011) with aCV4SEM:::BFGS approximated Hessian
#' @param SEM model of class Rcpp_SEMCpp. 
#' @param regularizedParameterLabels labels of regularized parameters
#' @param lambda lasso tuning parameter. Higher values = higher penalty
#' @param alpha 0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge, in between is elastic net
#' @param adaptiveLassoWeights vector with weights for adaptive LASSO. Set to NULL if not using adaptive LASSO
#' @param initialHessian option to provide an initial Hessian to the optimizer
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
#' @param c1 c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param epsOut Stopping criterion for outer iterations
#' @param epsIn Stopping criterion for inner iterations
#' @param useMultipleConvergencCriteria if set to TRUE, GLMNET will also check the change in fit and the change in parameters. If any convergence criterion is met, the optimization stops
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
GLMNET <- function(SEM, 
                   regularizedParameterLabels, 
                   lambda, 
                   alpha,
                   adaptiveLassoWeights,
                   initialHessian,
                   stepSize,
                   c1,
                   c2,
                   sig,
                   gam,
                   maxIterOut,
                   maxIterIn,
                   maxIterLine,
                   epsOut,
                   epsIn,
                   useMultipleConvergencCriteria,
                   verbose = 0){
  if(!((0 <= alpha) && (alpha <= 1))) stop("alpha must be in [0,1]")
  
  # Setup
  N <- nrow(SEM$rawData)
  
  # save current state
  initialParameters <- aCV4SEM:::getParameters(SEM, raw = TRUE)
  if(is.null(initialHessian)){
    initialHessian <- aCV4SEM:::getHessian(SEM = SEM, raw = TRUE)
  }
  initialGradients <- try(aCV4SEM:::getGradients(SEM = SEM, raw = TRUE))
  if(is(initialGradients, "try-error")){stop("Could not compute gradients at initial parameter values.")}
  
  # initialize parameter values
  newParameters <- initialParameters
  parameterLabels <- names(newParameters)
  regularized <- parameterLabels%in%regularizedParameterLabels
  if(is.null(adaptiveLassoWeights)){
    if(alpha==1) {
      # lasso or elastic net
      adaptiveLassoWeights <- rep(1, length(newParameters))
    }else if(alpha==0){
      # ridge
      adaptiveLassoWeights <- rep(2, length(newParameters))
    }else{
      stop("adaptiveLassoWeights missing without default.")
    }
    names(adaptiveLassoWeights) <- parameterLabels
  }
  
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
    warning("Initial Hessian is not positive definite. Using diagonal matrix instead.")
    newHessian <- diag(100, nrow(initialHessian))
    rownames(newHessian) <- parameterLabels
    colnames(newHessian) <- parameterLabels
  }
  
  if(alpha != 1){
    # elastic net or ridge are used -> we add the differentiable part of 
    # the penalty to the gradients and the Hessian
    if(length(unique(adaptiveLassoWeights))!= 1) stop("ridge and elastic net can not be combined with adaptive lasso weights")
    # we multiply with the adaptive lasso weights. These are all 1 in case of elastic net. In case of ridge, they are
    # set to 0. This allows for one single unified optimization procedure 
    penaltyFunctionArguments <- list("regularizedParameterLabels" = regularizedParameterLabels,
                                     "lambda" = unique(adaptiveLassoWeights)*.5*lambda*(1-alpha)*N)
    newGradients <- newGradients + aCV4SEM::ridgeGradient(parameters = newParameters, 
                                                          penaltyFunctionArguments = penaltyFunctionArguments)[names(newGradients)]
    newHessian <- newHessian + ridgeHessian(parameters = newParameters, 
                                            penaltyFunctionArguments = penaltyFunctionArguments)[rownames(newHessian), colnames(newHessian)]
  }
  
  newM2LL <- SEM$m2LL
  newRegM2LL <- newM2LL + aCV4SEM::getPenaltyValue(parameters = newParameters,
                                                   lambda = N*lambda,
                                                   alpha = alpha,
                                                   regularizedParameterLabels = regularizedParameterLabels,
                                                   adaptiveLassoWeights = adaptiveLassoWeights)
  
  while(iterOut < maxIterOut){
    iterOut <- iterOut +1
    
    if(iterOut > 1){
      parameterChange <- (newParameters - oldParameters)/(newParameters+1e-9) # add small jitter as zeroed parameters will otherwise 
      # result in division by 0
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
      convergenceCriterion <- try(max(diag(diag(oldHessian))%*%direction^2) < epsOut, silent = TRUE)
      if(is(convergenceCriterion, "try-error") || is.na(convergenceCriterion)){
        converged <- FALSE
        # save last working parameters
        newParameters <- oldParameters
        newGradients <- oldGradients
        newHessian <- oldHessian
        warning(paste("The model did NOT CONVERGE for lambda = ", lambda, "."))
        break
      }
      if(convergenceCriterion){
        break
      }
      if(useMultipleConvergencCriteria && abs(regM2LLChange) < 1e-8 && max(abs(parameterChange)) < 1e-8){
        # second convergence criterion: no more change in fit
        break
      }
    }
    
    
    ## inner loop: optimize directions
    direction <- aCV4SEM:::innerGLMNET(parameters = oldParameters,
                                       subGroupGradient = oldGradients, 
                                       subGroupHessian = oldHessian, 
                                       subGroupLambda = N*lambda*alpha, # update for lasso part
                                       regularized = regularized, 
                                       adaptiveLassoWeights = adaptiveLassoWeights,
                                       maxIter = maxIterIn, 
                                       epsBreak = epsIn
    )
    direction <- c(direction)
    names(direction) <- parameterLabels
    
    # perform Line Search
    
    stepSize_k <- aCV4SEM:::GLMNETLineSearch(SEM = SEM, 
                                             adaptiveLassoWeights = adaptiveLassoWeights, 
                                             parameterLabels = parameterLabels,
                                             regularizedParameterLabels = regularizedParameterLabels,
                                             lambda = N*lambda,
                                             alpha = alpha,
                                             oldParameters = oldParameters, 
                                             oldM2LL = oldM2LL, 
                                             oldGradients = oldGradients,
                                             oldHessian = oldHessian,
                                             direction = direction,
                                             stepSize= stepSize, 
                                             sig = sig,
                                             gam = gam, 
                                             maxIterLine = maxIterLine)
    
    newParameters <- oldParameters+stepSize_k*direction
    
    # update model: set parameter values and compute
    SEM <- aCV4SEM:::setParameters(SEM = SEM, labels = names(newParameters), values = newParameters, raw = TRUE)
    SEM <- fit(SEM)
    
    # get fit
    newM2LL <- SEM$m2LL
    newRegM2LL <- newM2LL + aCV4SEM::getPenaltyValue(parameters = newParameters,
                                                     lambda = N*lambda,
                                                     alpha = alpha,
                                                     regularizedParameterLabels = regularizedParameterLabels,
                                                     adaptiveLassoWeights = adaptiveLassoWeights)
    
    # extract gradients:
    newGradients <- try(aCV4SEM:::getGradients(SEM = SEM, raw = TRUE))
    if(is(newGradients, "try-error")){
      stop("Could not compute gradients at new location.")
    }
    
    if(alpha != 1){
      # add derivative of differentiable part of the penalty:
      newGradients <- newGradients + aCV4SEM::ridgeGradient(parameters = newParameters, 
                                                            penaltyFunctionArguments = penaltyFunctionArguments)[names(newGradients)]
    }
    
    # Approximate Hessian using aCV4SEM:::BFGS
    newHessian <- aCV4SEM:::BFGS(oldParameters = oldParameters, 
                                 oldGradients = oldGradients, 
                                 oldHessian = oldHessian, 
                                 newParameters = newParameters, 
                                 newGradients = newGradients)
    
    if(verbose > 0){
      
      cat(paste0("\r",
                 "## [", iterOut,
                 "] m2LL: ", sprintf('%.3f',newM2LL),
                 " | regM2LL:  ", sprintf('%.3f',newRegM2LL),
                 " | zeroed: ", sum(newParameters[regularizedParameterLabels] == 0),
                 " ##")
      )
      flush.console()
    }
    
  }
  # warnings
  if(iterOut == maxIterOut){
    warning(paste("For lambda = ", lambda, "the maximum number of iterations was reached. Try with a higher maxIterOut or with smaller lambda-steps."))
  }
  
  if(alpha != 1){
    # remove derivative of differentiable part of the penalty; we want the gradients of the log-Likelihood, not those including the penalty
    newGradients <- newGradients - aCV4SEM::ridgeGradient(parameters = newParameters, 
                                                          penaltyFunctionArguments = penaltyFunctionArguments)[names(newGradients)]
    # remove hessian of differentiable part of the penalty; we want the hessian of the log-Likelihood, not the one including the penalty
    newHessian <- newHessian - aCV4SEM::ridgeHessian(parameters = newParameters, 
                                                     penaltyFunctionArguments = penaltyFunctionArguments)[rownames(newHessian), colnames(newHessian)]
  }
  
  return(list("SEM" = SEM, 
              "parameters" = newParameters, 
              "m2LL" = newM2LL,
              "regM2LL" = newRegM2LL,
              "nonZeroParameters" = sum(newParameters != 0),
              "gradients" = newGradients, 
              "Hessian" = newHessian, 
              "convergence" = converged))
}




#### Line search ####

#' GLMNETLineSearch
#'
#' performs the line search procedure described by Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421 Equation 20.
#'
#'
#' @param SEM model of class Rcpp_SEMCpp. 
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
GLMNETLineSearch <- function(SEM, 
                             regularizedParameterLabels,
                             adaptiveLassoWeights, 
                             parameterLabels,
                             lambda,
                             alpha,
                             oldParameters, 
                             oldM2LL, 
                             oldGradients,
                             oldHessian,
                             direction,
                             stepSize, 
                             sig,
                             gam, 
                             maxIterLine){
  
  # get penalized M2LL for step size 0:
  pen_0 <- aCV4SEM::getPenaltyValue(parameters = oldParameters,
                                    lambda = lambda,
                                    alpha = alpha,
                                    regularizedParameterLabels = regularizedParameterLabels,
                                    adaptiveLassoWeights = adaptiveLassoWeights)
  f_0 <- oldM2LL + pen_0
  
  i <- 0
  stepSizeInit <- stepSize
  if(stepSizeInit >= 1) stepSizeInit <- .9
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
    p_new <- aCV4SEM::getPenaltyValue(parameters = newParameters,
                                      lambda = lambda,
                                      alpha = alpha,
                                      regularizedParameterLabels = regularizedParameterLabels,
                                      adaptiveLassoWeights = adaptiveLassoWeights)
    f_new <- newM2LL + p_new
    
    # test line search criterion
    lineCriterion <- f_new - f_0 <= sig*stepSize*(t(oldGradients)%*%direction + gam*t(direction)%*%oldHessian%*%direction + p_new - pen_0)
    if(lineCriterion){
      break
    }
    i <- i+1
    if(i >= maxIterLine){
      break
    }
  }
  return(stepSize)
}

#' BFGS
#'
#' computes the BFGS Hessian approximation
#'
#'
#' @param oldParameters parameters of previous iteration
#' @param oldGradients gradients of previous iteration
#' @param oldHessian Hessian of previous iteration
#' @param newParameters parameters of current iteration
#' @param newGradients gradients of current iteration
#' @param cautious boolean: should the update be skipped if it would result in a non positive definite Hessian?
#' @param hessianEps controls when the update of the Hessian approximation is skipped
BFGS <- function(oldParameters, 
                 oldGradients, 
                 oldHessian, 
                 newParameters, 
                 newGradients, 
                 cautious = TRUE, 
                 hessianEps = .001){
  
  y <- newGradients-oldGradients
  d <- newParameters-oldParameters
  
  # test if positive definiteness is ensured
  skipUpdate <- try((t(y)%*%d < hessianEps) && cautious, silent = TRUE)
  if(any(class(skipUpdate) == "try-error") || skipUpdate || is.na(skipUpdate)){
    # Hessian might become non-positive definite. Return without update
    return(oldHessian)
  }
  if(t(y)%*%d < 0){
    warning("Hessian update possibly non-positive definite.")
  }
  
  newHessian <- oldHessian - (oldHessian%*%d%*%t(d)%*%oldHessian)/as.numeric(t(d)%*%oldHessian%*%d) + (y%*%t(y))/as.numeric(t(y)%*%d)
  
  if(anyNA(newHessian)){
    warning("Invalid Hessian. Returning previous Hessian")
    return(oldHessian)
  }
  HessianEigen <- eigen(newHessian)
  iscomplex <- any(!Im(HessianEigen$values) == 0)
  if(iscomplex){
    eigenVectors <- Re(HessianEigen$vectors)
    eigenValues <- Re(HessianEigen$values)
    newHessian <- eigenVectors%*%diag(eigenValues)%*%t(eigenVectors)
  }
  if(any(HessianEigen$values < 0)){
    while(any(eigen(newHessian)$values < 0)){
      newHessian <- newHessian + .01*diag(nrow(newHessian))
    }
  }
  
  return(newHessian)
  
}

#' getPenaltyValue
#'
#' returns the penalty value at a specific location
#'
#'
#' @param parameters vector with labeled parameter values
#' @param lambda lasso tuning parameter. Higher values = higher penalty
#' @param alpha 0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge
#' @param regularizedParameterLabels labels of regularized parameters
#' @param adaptiveLassoWeights vector with weights for adaptive LASSO. Set to vector of ones if not using adaptive LASSO
#' @export
getPenaltyValue <- function(parameters,
                            lambda,
                            alpha,
                            regularizedParameterLabels,
                            adaptiveLassoWeights){
  return(alpha*lambda*sum(adaptiveLassoWeights[regularizedParameterLabels]*abs(parameters[regularizedParameterLabels])) + # lasso
           .5*(1-alpha)*lambda*sum(adaptiveLassoWeights[regularizedParameterLabels]*parameters[regularizedParameterLabels]^2) # ridge
  )
}
