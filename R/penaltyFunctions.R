#### RIDGE ####

ridge <- function(parameters, 
                  penaltyFunctionArguments){
  return(penaltyFunctionArguments$lambda*sum(parameters[penaltyFunctionArguments$regularizedParameterLabels]^2))
}

ridgeGradient <- function(parameters, 
                          penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- penaltyFunctionArguments$lambda*2*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]
  return(gradients)
  
}

ridgeHessian <- function(parameters, 
                         penaltyFunctionArguments){
  hessian <- matrix(0, 
                    length(parameters),
                    length(parameters),
                    dimnames = list(names(parameters),
                                    names(parameters)))
  diag(hessian[penaltyFunctionArguments$regularizedParameterLabels,
               penaltyFunctionArguments$regularizedParameterLabels]) <- 2*penaltyFunctionArguments$lambda
  return(hessian)
}

#### LASSO ####

smoothLASSO <- function(parameters, 
                        penaltyFunctionArguments
){
  return(penaltyFunctionArguments$lambda*sum(sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels ])^2 + penaltyFunctionArguments$eps)))
}

smoothLASSOGradient <- function(parameters, 
                                penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- penaltyFunctionArguments$lambda*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(gradients)
}

smoothLASSOHessian <- function(parameters, 
                               penaltyFunctionArguments){
  hessian <- matrix(0, 
                    length(parameters),
                    length(parameters),
                    dimnames = list(names(parameters),
                                    names(parameters)))
  diag(hessian[penaltyFunctionArguments$regularizedParameterLabels,
               penaltyFunctionArguments$regularizedParameterLabels]) <- -penaltyFunctionArguments$lambda*
    (parameters[penaltyFunctionArguments$regularizedParameterLabels])^2*
    ((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps)^(-3/2) +
    penaltyFunctionArguments$lambda*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(hessian)
}

#### ADAPTIVE LASSO ####
smoothAdaptiveLASSO <- function(parameters, 
                                penaltyFunctionArguments){
  return(sum(penaltyFunctionArguments$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
               sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps)))
}

smoothAdaptiveLASSOGradient <- function(parameters, 
                                        penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- penaltyFunctionArguments$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(gradients)
}

smoothAdaptiveLASSOHessian <- function(parameters, 
                                       penaltyFunctionArguments){
  hessian <- matrix(0, 
                    length(parameters),
                    length(parameters),
                    dimnames = list(names(parameters),
                                    names(parameters)))
  diag(hessian[penaltyFunctionArguments$regularizedParameterLabels,
               penaltyFunctionArguments$regularizedParameterLabels]) <- -penaltyFunctionArguments$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
    (parameters[penaltyFunctionArguments$regularizedParameterLabels])^2*
    ((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps)^(-3/2) +
    penaltyFunctionArguments$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(hessian)
}

#### ELASTIC NET ####

smoothElasticNet <- function(parameters, 
                             penaltyFunctionArguments){
  # using the implementation from https://github.com/psyphh/lslx/blob/ver.0.6.11/src/lslxOptimizer.cpp
  lassoPart <- smoothLASSO(parameters = parameters,
                           penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                           "lambda" = penaltyFunctionArguments$alpha*penaltyFunctionArguments$lambda,
                                                           "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- ridge(parameters = parameters,
                     penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                     "lambda" = 0.5*penaltyFunctionArguments$lambda*(1-penaltyFunctionArguments$alpha)
                     )
  )
  return(lassoPart + ridgePart)
}

smoothElasticNetGradient <- function(parameters, 
                                     penaltyFunctionArguments){
  # using the implementation from https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
  lassoPart <- smoothLASSOGradient(parameters = parameters,
                                   penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                                   "lambda" = penaltyFunctionArguments$alpha*penaltyFunctionArguments$lambda,
                                                                   "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- ridgeGradient(parameters = parameters,
                             penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                             "lambda" = 0.5*penaltyFunctionArguments$lambda*(1-penaltyFunctionArguments$alpha)
                             )
  )
  return(lassoPart + ridgePart)
}

smoothElasticNetHessian <- function(parameters, 
                                     penaltyFunctionArguments){
  # using the implementation from https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
  lassoPart <- smoothLASSOHessian(parameters = parameters,
                                   penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                                   "lambda" = penaltyFunctionArguments$alpha*penaltyFunctionArguments$lambda,
                                                                   "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- ridgeHessian(parameters = parameters,
                             penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                             "lambda" = 0.5*penaltyFunctionArguments$lambda*(1-penaltyFunctionArguments$alpha)
                             )
  )
  return(lassoPart + ridgePart)
}

#### MCP ####

smoothMCP <- function(parameters, 
                      regularizedParameterLabels, 
                      lambda, 
                      eps = 1e-4){
  stop("Not yet implmented")
}

#### SCAD ####
smoothSCAD <- function(parameters, 
                       regularizedParameterLabels, 
                       lambda, 
                       eps = 1e-4){
  stop("Not yet implmented")
}