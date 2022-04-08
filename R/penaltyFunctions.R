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

ridge <- function(parameters, 
                  regularizedParameterLabels, 
                  lambda){
  return(lambda*sum(parameters[regularizedParameterLabels]^2))
}

smoothAdaptiveLASSO <- function(parameters, 
                                regularizedParameterLabels,
                                lambdas, 
                                eps = 1e-4){
  return(sum(lambdas[regularizedParameterLabels]*sqrt((parameters[regularizedParameterLabels])^2 + eps)))
}


smoothElasticNet <- function(parameters, 
                             regularizedParameterLabels, 
                             lambda, 
                             alpha,
                             eps = 1e-4){
  # using the implementation from https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
  return(alpha*lambda*sum(sqrt((parameters[regularizedParameterLabels])^2 + eps))
         + 0.5*alpha*(1-lambda)*sum(parameters[regularizedParameterLabels]^2)
  )
}

smoothMCP <- function(parameters, 
                      regularizedParameterLabels, 
                      lambda, 
                      eps = 1e-4){
  stop("Not yet implmented")
}

smoothSCAD <- function(parameters, 
                       regularizedParameterLabels, 
                       lambda, 
                       eps = 1e-4){
  stop("Not yet implmented")
}