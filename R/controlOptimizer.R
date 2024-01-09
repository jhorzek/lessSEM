#' controlIsta
#' 
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param saveDetails when set to TRUE, additional details about the individual
#' models are save. Currently, this are the implied means and covariances.
#'  Note: This may take a lot of memory!
#' @param L0 L0 controls the step size used in the first iteration
#' @param eta eta controls by how much the step size changes in the
#' inner iterations with (eta^i)*L, where i is the inner iteration
#' @param accelerate boolean: Should the acceleration outlined in 
#' Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and Trends 
#' in Optimization, 1(3), 123–231., p. 152 be used?
#' @param maxIterOut maximal number of outer iterations
#' @param maxIterIn maximal number of inner iterations
#' @param breakOuter change in fit required to break the outer iteration. Note: The
#' value will be multiplied internally with sample size N as the -2log-Likelihood
#' depends directly on the sample size
#' @param convCritInner this is related to the inner breaking condition.
#' 0 = ista, as presented by Beck & Teboulle (2009); see Remark 3.1 on p. 191 (ISTA with backtracking)
#' 1 = gist, as presented by Gong et al. (2013) (Equation 3)
#' @param sigma sigma in (0,1) is used by the gist convergence criterion. larger
#'  sigma enforce larger improvement in fit
#' @param stepSizeInheritance how should step sizes be carried forward from iteration to iteration? 
#' 0 = resets the step size to L0 in each iteration
#' 1 = takes the previous step size as initial value for the next iteration
#' 3 = Barzilai-Borwein procedure 
#' 4 = Barzilai-Borwein procedure, but sometimes resets the step size; this can help when the optimizer is caught in a bad spot.
#' @param verbose if set to a value > 0, the fit every "verbose" iterations is printed.
#' @param nCores number of core to use. Multi-core support is provided by RcppParallel and only supported for SEM, not for general purpose optimization.
#' @returns object of class controlIsta
#' @examples
#' control <- controlIsta()
#' @export
controlIsta <- function(
    startingValues = "est",
    saveDetails = FALSE,
    L0 = .1,
    eta = 2,
    accelerate = TRUE,
    maxIterOut = 10000,
    maxIterIn = 1000,
    breakOuter = .00000001,
    convCritInner = 1,
    sigma = .1,
    stepSizeInheritance = ifelse(accelerate,1,3),
    verbose = 0,
    nCores = 1
){
  control <- as.list(environment())
  class(control) <- "controlIsta"
  return(control)
}

#' controlGlmnet
#' 
#' Control the GLMNET optimizer.
#' 
#' @param startingValues option to provide initial starting values. Only used 
#' for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values 
#' of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param initialHessian option to provide an initial Hessian to the optimizer. 
#' Must have row and column names corresponding to the parameter labels. use 
#' getLavaanParameters(lavaanModel) to 
#' see those labels. If set to "gradNorm", the maximum of the gradients at the starting
#' values times the stepSize will be used. This is adapted from Optim.jl
#' https://github.com/JuliaNLSolvers/Optim.jl/blob/f43e6084aacf2dabb2b142952acd3fbb0e268439/src/multivariate/solvers/first_order/bfgs.jl#L104
#' If set to "compute", the initial hessian will be computed. If set to a single 
#' value, a diagonal matrix with the single value along the diagonal will be used.
#' The default is "lavaan" which extracts the Hessian from the lavaanModel. This Hessian
#' will typically deviate from that of the internal SEM represenation of lessSEM (due to
#' the transformation of the variances), but works quite well in practice.
#' @param saveDetails when set to TRUE, additional details about the individual
#' models are save. Currently, this are the Hessian and the implied means and covariances.
#'  Note: This may take a lot of memory!
#' @param stepSize Initial stepSize of the outer iteration 
#' (theta_next = theta_previous + stepSize * Stepdirection)
#' @param sigma only relevant when lineSearch = 'GLMNET'. Controls the sigma 
#' parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET 
#' for l1-regularized logistic regression. The Journal of Machine Learning Research, 
#' 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gamma Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
#' An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning 
#' Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param breakOuter Stopping criterion for outer iterations
#' @param breakInner Stopping criterion for inner iterations
#' @param convergenceCriterion which convergence criterion should be used for the outer iterations? possible are 0 = GLMNET, 1 = fitChange, 2 = gradients.
#' Note that in case of gradients and GLMNET, we divide the gradients (and the Hessian) of the log-Likelihood by N as it would otherwise be
#' considerably more difficult for larger sample sizes to reach the convergence criteria.
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
#' @param nCores number of core to use. Multi-core support is provided by RcppParallel and only supported for SEM, not for general purpose optimization.
#' @returns object of class controlGlmnet
#' @examples
#' control <- controlGlmnet()
#' @export
controlGlmnet <- function(
    startingValues = "est",
    initialHessian = ifelse(all(startingValues=="est"),"lavaan","compute"),
    saveDetails = FALSE,
    stepSize = .9,
    sigma = 1e-5,
    gamma = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 500,
    breakOuter = 1e-8,
    breakInner = 1e-10,
    convergenceCriterion = 0,
    verbose = 0,
    nCores = 1
){
  control <- as.list(environment())
  class(control) <- "controlGlmnet"
  return(control)
}

#' controlBFGS
#' 
#' Control the BFGS optimizer.
#' 
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param initialHessian option to provide an initial Hessian to the optimizer. Must have row and column names corresponding to the parameter labels. use getLavaanParameters(lavaanModel) to 
#' see those labels. If set to "gradNorm", the maximum of the gradients at the starting
#' values times the stepSize will be used. This is adapted from Optim.jl
#' https://github.com/JuliaNLSolvers/Optim.jl/blob/f43e6084aacf2dabb2b142952acd3fbb0e268439/src/multivariate/solvers/first_order/bfgs.jl#L104
#' If set to a single value, a diagonal matrix with the single value along the diagonal will be used.
#' The default is "lavaan" which extracts the Hessian from the lavaanModel. This Hessian
#' will typically deviate from that of the internal SEM represenation of lessSEM (due to
#' the transformation of the variances), but works quite well in practice.
#' @param saveDetails when set to TRUE, additional details about the individual
#' models are save. Currently, this are the Hessian and the implied means and covariances.
#'  Note: This may take a lot of memory!
#' @param stepSize Initial stepSize of the outer iteration (theta_next = theta_previous + stepSize * Stepdirection)
#' @param sigma only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gamma Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param breakOuter Stopping criterion for outer iterations
#' @param breakInner Stopping criterion for inner iterations
#' @param convergenceCriterion which convergence criterion should be used for the outer iterations? possible are 0 = GLMNET, 1 = fitChange, 2 = gradients.
#' Note that in case of gradients and GLMNET, we divide the gradients (and the Hessian) of the log-Likelihood by N as it would otherwise be
#' considerably more difficult for larger sample sizes to reach the convergence criteria.
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
#' @param nCores number of core to use. Multi-core support is provided by RcppParallel and only supported for SEM, not for general purpose optimization.
#' @returns object of class controlBFGS
#' @examples
#' control <- controlBFGS()
#' @export
controlBFGS <- function(
    startingValues = "est",
    initialHessian = ifelse(all(startingValues=="est"),"lavaan","compute"),
    saveDetails = FALSE,
    stepSize = .9,
    sigma = 1e-5,
    gamma = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 500,
    breakOuter = 1e-8,
    breakInner = 1e-10,
    convergenceCriterion = 0,
    verbose = 0,
    nCores = 1
){
  control <- as.list(environment())
  class(control) <- "controlBFGS"
  return(control)
}

#' .setupMulticore
#' 
#' setup for multi-core support
#' @param control object created with controlBFGS, controlIsta or controlGlmnet function
#' @return nothing
#' @keywords internal
.setupMulticore <- function(control){
  if(RcppParallel::defaultNumThreads() < control$nCores)
    warning("Your selected number of cores (", control$nCores,
            ") is larger than the number of cores detected by RcppParallel (",
            RcppParallel::defaultNumThreads(), "). You may consider using fewer cores."
    )
  
  RcppParallel::setThreadOptions(numThreads = control$nCores)
  
}


#' .adaptBreakingForWls
#' 
#' wls needs smaller breaking points than ml
#' @param lavaanModel single model or vector of models
#' @param currentBreaking current breaking condition value
#' @param selectedDefault was default breaking condition selected?
#' @return updated breaking
.adaptBreakingForWls <- function(lavaanModel, currentBreaking, selectedDefault){
  
  if(is.vector(lavaanModel)){
    for(i in 1:length(lavaanModel)){
      currentBreaking <- min(currentBreaking, 
                             .adaptBreakingForWls(lavaanModel = lavaanModel[[i]], 
                                                  currentBreaking = currentBreaking,
                                                  selectedDefault = selectedDefault))
    }
    return(currentBreaking)
  }
  
  if(selectedDefault && (tolower(lavaanModel@Options$estimator) %in% c("uls","wls", "dwls", "gls"))){
    return(1e-11)
  }else{
    return(currentBreaking)
  }
}