#' controlGLMNET
#' 
#' Control the GLMNET optimizer.
#' 
#' @param startingValues option to provide initial starting values. Only used for the first lambda.
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
#' @export
controlGLMNET <- function(
  startingValues = NULL,
  initialHessian = NULL,
  stepSize = 1,
  c1 = 1e-04,
  c2 = 0.9,
  sig = 10^(-5),
  gam = 0,
  maxIterOut = 100,
  maxIterIn = 1000,
  maxIterLine = 500,
  epsOut = 1e-06,
  epsIn = 1e-06,
  useMultipleConvergencCriteria = TRUE,
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
    useMultipleConvergencCriteria = useMultipleConvergencCriteria,
    verbose = verbose
    )
  )
}