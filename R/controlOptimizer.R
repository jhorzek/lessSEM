#' controlGLMNET
#' 
#' Control the GLMNET optimizer.
#' 
#' @param addMeans If lavaanModel has meanstructure = FALSE, addMeans = TRUE will add a mean structure. FALSE will set the means of the observed variables to the average
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param initialHessian option to provide an initial Hessian to the optimizer. Must have row and column names corresponding to the parameter labels. use getLavaanParameters(lavaanModel) to 
#' see those labels. If set to "scoreBased", the outer product of the scores will be used as an approximation (see https://en.wikipedia.org/wiki/Berndt%E2%80%93Hall%E2%80%93Hall%E2%80%93Hausman_algorithm).
#' If set to "compute", the initial hessian will be computed. If set to a single value, a diagonal matrix with the single value along the diagonal will be used.
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param epsOut Stopping criterion for outer iterations
#' @param epsIn Stopping criterion for inner iterations
#' @param convergenceCriterion which convergence criterion should be used for the outer iterations? possible are "GLMNET", "gradients", "fitChange"
#' @param saveHessian should the Hessian be saved for later use? Note: This may take a lot of memory!
#' @param activeSet Option to only use a subset of the individuals in the data set. Logical vector of length N indicating which subjects should remain in the sample.
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
#' @export
controlGLMNET <- function(
    addMeans = TRUE,
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
    epsOut = 1e-8,
    epsIn = 1e-10,
    convergenceCriterion = "GLMNET",
    saveHessian = FALSE,
    verbose = 0
){
  return(
    list(
      addMeans = addMeans,
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
      convergenceCriterion = convergenceCriterion,
      saveHessian = saveHessian,
      verbose = verbose
    )
  )
}

#' controlQuasiNewtonBFGS
#' 
#' Control the quasiNewtonBFGS optimizer.
#' 
#' @param addMeans If lavaanModel has meanstructure = FALSE, addMeans = TRUE will add a mean structure. FALSE will set the means of the observed variables to the average
#' @param startingValues option to provide initial starting values. Only used for the first lambda. Three options are supported. Setting to "est" will use the estimates
#' from the lavaan model object. Setting to "start" will use the starting values of the lavaan model. Finally, a labeled vector with parameter
#' values can be passed to the function which will then be used as starting values.
#' @param initialHessian option to provide an initial Hessian to the optimizer. Must have row and column names corresponding to the parameter labels. use getLavaanParameters(lavaanModel) to 
#' see those labels. If set to "scoreBased", the outer product of the scores will be used as an approximation (see https://en.wikipedia.org/wiki/Berndt%E2%80%93Hall%E2%80%93Hall%E2%80%93Hausman_algorithm).
#' If set to "compute", the initial hessian will be computed. If set to a single value, a diagonal matrix with the single value along the diagonal will be used.
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIterOut Maximal number of outer iterations
#' @param maxIterIn Maximal number of inner iterations
#' @param maxIterLine Maximal number of iterations for the line search procedure
#' @param epsOut Stopping criterion for outer iterations
#' @param epsIn Stopping criterion for inner iterations
#' @param convergenceCriterion which convergence criterion should be used for the outer iterations? possible are "GLMNET", "gradients", "fitChange"
#' @param saveHessian should the Hessian be saved for later use? Note: This may take a lot of memory!
#' @param activeSet Option to only use a subset of the individuals in the data set. Logical vector of length N indicating which subjects should remain in the sample.
#' @param verbose 0 prints no additional information, > 0 prints GLMNET iterations
#' @export
controlQuasiNewtonBFGS <- function(
    addMeans = TRUE,
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
    epsOut = 1e-8,
    epsIn = 1e-10,
    convergenceCriterion = "GLMNET",
    saveHessian = FALSE,
    verbose = 0
){
  return(
    list(
      addMeans = addMeans,
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
      convergenceCriterion = convergenceCriterion,
      saveHessian = saveHessian,
      verbose = verbose
    )
  )
}