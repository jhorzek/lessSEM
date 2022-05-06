#### RIDGE ####

#' ridge
#' 
#' ridge penalty function
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value) 
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters)
#' @export
ridge <- function(parameters, 
                  tuningParameters,
                  penaltyFunctionArguments){
  return(tuningParameters$lambda*sum(parameters[penaltyFunctionArguments$regularizedParameterLabels]^2))
}

#' ridgeGradient
#' 
#' ridge gradient function
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value) 
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters)
#' @export
ridgeGradient <- function(parameters, 
                          tuningParameters,
                          penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- tuningParameters$lambda*2*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]
  return(gradients)
  
}

#' ridgeGradient
#' 
#' ridge Hessian function
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value) 
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters)
#' @export
ridgeHessian <- function(parameters, 
                         tuningParameters,
                         penaltyFunctionArguments){
  hessian <- matrix(0, 
                    length(parameters),
                    length(parameters),
                    dimnames = list(names(parameters),
                                    names(parameters)))
  diag(hessian[penaltyFunctionArguments$regularizedParameterLabels,
               penaltyFunctionArguments$regularizedParameterLabels]) <- 2*tuningParameters$lambda
  return(hessian)
}

#### LASSO ####
#' smoothLASSO
#' 
#' smoothed version of non-differentiable LASSO penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothLASSO <- function(parameters, 
                        tuningParameters,
                        penaltyFunctionArguments
){
  return(tuningParameters$lambda*sum(sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels ])^2 + penaltyFunctionArguments$eps)))
}

#' smoothLASSOGradient
#' 
#' smoothed version of non-differentiable LASSO gradient
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothLASSOGradient <- function(parameters, 
                                tuningParameters,
                                penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- tuningParameters$lambda*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(gradients)
}

#' smoothLASSOHessian
#' 
#' smoothed version of non-differentiable LASSO Hessian
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothLASSOHessian <- function(parameters,
                               tuningParameters,
                               penaltyFunctionArguments){
  hessian <- matrix(0, 
                    length(parameters),
                    length(parameters),
                    dimnames = list(names(parameters),
                                    names(parameters)))
  diag(hessian[penaltyFunctionArguments$regularizedParameterLabels,
               penaltyFunctionArguments$regularizedParameterLabels]) <- -tuningParameters$lambda*
    (parameters[penaltyFunctionArguments$regularizedParameterLabels])^2*
    ((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps)^(-3/2) +
    tuningParameters$lambda*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(hessian)
}

#### ADAPTIVE LASSO ####
#
#' smoothAdaptiveLASSO
#' 
#' smoothed version of non-differentiable adaptive LASSO penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambdas (vector with one tuning parameter value for each parameter)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothAdaptiveLASSO <- function(parameters, 
                                tuningParameters,
                                penaltyFunctionArguments){
  return(sum(tuningParameters$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
               sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps)))
}

#' smoothAdaptiveLASSOGradient
#' 
#' smoothed version of non-differentiable adaptive LASSO gradient
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambdas (vector with one tuning parameter value for each parameter)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothAdaptiveLASSOGradient <- function(parameters, 
                                        tuningParameters,
                                        penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- tuningParameters$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(gradients)
}

#' smoothAdaptiveLASSOHessian
#' 
#' smoothed version of non-differentiable adaptive LASSO Hessian
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambdas (vector with one tuning parameter value for each parameter)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothAdaptiveLASSOHessian <- function(parameters, 
                                       tuningParameters,
                                       penaltyFunctionArguments){
  hessian <- matrix(0, 
                    length(parameters),
                    length(parameters),
                    dimnames = list(names(parameters),
                                    names(parameters)))
  diag(hessian[penaltyFunctionArguments$regularizedParameterLabels,
               penaltyFunctionArguments$regularizedParameterLabels]) <- -tuningParameters$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
    (parameters[penaltyFunctionArguments$regularizedParameterLabels])^2*
    ((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps)^(-3/2) +
    tuningParameters$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(hessian)
}

#### ELASTIC NET ####
#' smoothElasticNet
#' 
#' smoothed version of non-differentiable elastic LASSO penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambda (tuning parameter value), alpha (0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothElasticNet <- function(parameters, 
                             tuningParameters,
                             penaltyFunctionArguments){
  # using the implementation from https://github.com/psyphh/lslx/blob/ver.0.6.11/src/lslxOptimizer.cpp
  lassoPart <- smoothLASSO(parameters = parameters,
                           penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                           "lambda" = tuningParameters$alpha*tuningParameters$lambda,
                                                           "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- ridge(parameters = parameters,
                     penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                     "lambda" = tuningParameters$lambda*(1-tuningParameters$alpha)
                     )
  )
  return(lassoPart + ridgePart)
}

#' smoothElasticNetGradient
#' 
#' smoothed version of non-differentiable elastic LASSO gradient
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambda (tuning parameter value), alpha (0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothElasticNetGradient <- function(parameters, 
                                     tuningParameters,
                                     penaltyFunctionArguments){
  # using the implementation from https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
  lassoPart <- smoothLASSOGradient(parameters = parameters,
                                   penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                                   "lambda" = tuningParameters$alpha*tuningParameters$lambda,
                                                                   "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- ridgeGradient(parameters = parameters,
                             penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                             "lambda" = tuningParameters$lambda*(1-tuningParameters$alpha)
                             )
  )
  return(lassoPart + ridgePart)
}

#' smoothElasticNetHessian
#' 
#' smoothed version of non-differentiable elastic LASSO Hessian
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambda (tuning parameter value), alpha (0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @export
smoothElasticNetHessian <- function(parameters, 
                                    tuningParameters,
                                    penaltyFunctionArguments){
  # using the implementation from https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
  lassoPart <- smoothLASSOHessian(parameters = parameters,
                                  penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                                  "lambda" = tuningParameters$alpha*tuningParameters$lambda,
                                                                  "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- ridgeHessian(parameters = parameters,
                            penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                            "lambda" = tuningParameters$lambda*(1-tuningParameters$alpha)
                            )
  )
  return(lassoPart + ridgePart)
}
