smoothLASSO <- function(parameters, 
                        regularizedParameterLabels, 
                        lambda, 
                        eps = 1e-4){
  return(lambda*sum(sqrt((parameters[regularizedParameterLabels])^2 + eps)))
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