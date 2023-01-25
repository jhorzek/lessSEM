#' controlIsta
#' 
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
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
#' @returns object of class controlIsta
#' @export
controlIsta <- function(
    startingValues = "est",
    L0 = .1,
    eta = 2,
    accelerate = TRUE,
    maxIterOut = 10000,
    maxIterIn = 1000,
    breakOuter = .00000001,
    convCritInner = 1,
    sigma = .1,
    stepSizeInheritance = ifelse(accelerate,1,3),
    verbose = 0
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
#' see those labels. If set to "scoreBased", the outer product of the scores 
#' will be used as an approximation 
#' (see https://en.wikipedia.org/wiki/Berndt%E2%80%93Hall%E2%80%93Hall%E2%80%93Hausman_algorithm).
#' If set to "compute", the initial hessian will be computed. If set to a single 
#' value, a diagonal matrix with the single value along the diagonal will be used.
#' @param saveHessian should the Hessian be saved for later use? Note: This may take a lot of memory!
#' @param stepSize Initial stepSize of the outer iteration 
#' (theta_{k+1} = theta_k + stepSize * Stepdirection)
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
#' @returns object of class controlGlmnet
#' @export
controlGlmnet <- function(
    startingValues = "est",
    initialHessian = ifelse(all(startingValues=="est"),"compute",1),
    saveHessian = FALSE,
    stepSize = .9,
    sigma = 1e-5,
    gamma = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 500,
    breakOuter = 1e-8,
    breakInner = 1e-10,
    convergenceCriterion = 0,
    verbose = 0
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
#' see those labels. If set to "scoreBased", the outer product of the scores will be used as an approximation (see https://en.wikipedia.org/wiki/Berndt%E2%80%93Hall%E2%80%93Hall%E2%80%93Hausman_algorithm).
#' If set to "compute", the initial hessian will be computed. If set to a single value, a diagonal matrix with the single value along the diagonal will be used.
#' @param saveHessian should the Hessian be saved for later use? Note: This may take a lot of memory!
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
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
#' @returns object of class controlBFGS
#' @export
controlBFGS <- function(
    startingValues = "est",
    initialHessian = ifelse(all(startingValues=="est"),"compute",1),
    saveHessian = FALSE,
    stepSize = .9,
    sigma = 1e-5,
    gamma = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 500,
    breakOuter = 1e-8,
    breakInner = 1e-10,
    convergenceCriterion = 0,
    verbose = 0
){
  control <- as.list(environment())
  class(control) <- "controlBFGS"
  return(control)
}