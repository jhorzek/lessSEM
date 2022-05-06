#' controlGLMNET
#' 
#' Control the GLMNET optimizer.
#' 
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param initialHessian option to provide an initial Hessian to the optimizer. Must have row and column names corresponding to the parameter labels. use getLavaanParameters(lavaanModel) to 
#' see those labels. If set to "compute", the initial hessian will be computed. If set to a single value, a diagonal matrix with the single value along the diagonal will be used.
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param epsOut Stopping criterion for outer iterations
#' @param epsIn Stopping criterion for inner iterations
#' @param useMultipleConvergenceCriteria if set to TRUE, GLMNET will also check the change in fit and the change in parameters. If any convergence criterion is met, the optimization stops
#' @param regM2LLChangeEps if useMultipleConvergenceCriteria: change in fit which results in convergence
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
#' @export
controlGLMNET <- function(
  startingValues = "est",
  initialHessian = ifelse(all(startingValues=="est"),"compute",1),
  stepSize = 1,
  c1 = 1e-4,
  c2 = 0.9,
  sig = 1e-5,
  gam = 0,
  maxIterOut = 1000,
  maxIterIn = 1000,
  maxIterLine = 500,
  epsOut = 1e-5,
  epsIn = 1e-5,
  useMultipleConvergenceCriteria = TRUE,
  regM2LLChangeEps = 1e-6,
  verbose = 0
){
  return(
    list(
    startingValues = startingValues,
    initialHessian = initialHessian,
    stepSize = stepSize,
    c1 = c1,
    c2 = c2,
    sig = sig,
    gam = gam,
    maxIterOut = maxIterOut,
    maxIterIn = maxIterIn,
    maxIterLine = maxIterLine,
    epsOut = epsOut,
    epsIn = epsIn,
    useMultipleConvergenceCriteria = useMultipleConvergenceCriteria,
    regM2LLChangeEps = regM2LLChangeEps,
    verbose = verbose
    )
  )
}

#' controlQuasiNewtonBFGS
#' 
#' Control the quasiNewtonBFGS optimizer.
#' 
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param initialHessian option to provide an initial Hessian to the optimizer. Must have row and column names corresponding to the parameter labels. use getLavaanParameters(lavaanModel) to 
#' see those labels. If set to "compute", the initial hessian will be computed. If set to a single value, a diagonal matrix with the single value along the diagonal will be used.
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param epsOut Stopping criterion for outer iterations
#' @param epsIn Stopping criterion for inner iterations
#' @param useMultipleConvergenceCriteria if set to TRUE, GLMNET will also check the change in fit and the change in parameters. If any convergence criterion is met, the optimization stops
#' @param regM2LLChangeEps if useMultipleConvergenceCriteria: change in fit which results in convergence
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
#' @export
controlQuasiNewtonBFGS <- function(
    startingValues = "est",
    initialHessian = ifelse(all(startingValues=="est"),"compute",1),
    stepSize = 1,
    c1 = 1e-4,
    c2 = 0.9,
    sig = 1e-5,
    gam = 0,
    maxIterOut = 1000,
    maxIterIn = 1000,
    maxIterLine = 500,
    epsOut = 1e-5,
    epsIn = 1e-5,
    useMultipleConvergenceCriteria = TRUE,
    regM2LLChangeEps = 1e-6,
    verbose = 0
){
  return(
    list(
      startingValues = startingValues,
      initialHessian = initialHessian,
      stepSize = stepSize,
      c1 = c1,
      c2 = c2,
      sig = sig,
      gam = gam,
      maxIterOut = maxIterOut,
      maxIterIn = maxIterIn,
      maxIterLine = maxIterLine,
      epsOut = epsOut,
      epsIn = epsIn,
      useMultipleConvergenceCriteria = useMultipleConvergenceCriteria,
      regM2LLChangeEps = regM2LLChangeEps,
      verbose = verbose
    )
  )
}