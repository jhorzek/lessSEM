# The following functions are mostly used in tests to check the non-smooth implementations

#### RIDGE ####

#' .ridgeValue
#' 
#' ridge penalty function
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value) 
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters)
#' @returns penalty function value
#' @keywords internal
.ridgeValue <- function(parameters, 
                  tuningParameters,
                  penaltyFunctionArguments){
  return(tuningParameters$lambda*sum(parameters[penaltyFunctionArguments$regularizedParameterLabels]^2))
}

#' .ridgeGradient
#' 
#' ridge gradient function
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value) 
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters)
#' @returns gradient values
#' @keywords internal
.ridgeGradient <- function(parameters, 
                          tuningParameters,
                          penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- tuningParameters$lambda*2*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]
  return(gradients)
  
}

#' .ridgeHessian
#' 
#' ridge Hessian function
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value) 
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters)
#' @returns Hessian matrix
#' @keywords internal
.ridgeHessian <- function(parameters, 
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
#' .smoothLASSOValue
#' 
#' smoothed version of non-differentiable LASSO penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns penalty function value
#' @keywords internal
.smoothLASSOValue <- function(parameters, 
                        tuningParameters,
                        penaltyFunctionArguments
){
  return(tuningParameters$lambda*sum(sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels ])^2 + penaltyFunctionArguments$eps)))
}

#' .smoothLASSOGradient
#' 
#' smoothed version of non-differentiable LASSO gradient
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns gradient values
#' @keywords internal
.smoothLASSOGradient <- function(parameters, 
                                tuningParameters,
                                penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- tuningParameters$lambda*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(gradients)
}

#' .smoothLASSOHessian
#' 
#' smoothed version of non-differentiable LASSO Hessian
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns Hessian matrix
#' @keywords internal
.smoothLASSOHessian <- function(parameters,
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
#' .smoothAdaptiveLASSOValue
#' 
#' smoothed version of non-differentiable adaptive LASSO penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambdas (vector with one tuning parameter value for each parameter)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns penalty function value
#' @keywords internal
.smoothAdaptiveLASSOValue <- function(parameters, 
                                tuningParameters,
                                penaltyFunctionArguments){
  return(sum(tuningParameters$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
               sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps)))
}

#' .smoothAdaptiveLASSOGradient
#' 
#' smoothed version of non-differentiable adaptive LASSO gradient
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambdas (vector with one tuning parameter value for each parameter)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns gradient values
#' @keywords internal
.smoothAdaptiveLASSOGradient <- function(parameters, 
                                        tuningParameters,
                                        penaltyFunctionArguments){
  gradients <- rep(0, length(parameters))
  names(gradients) <- names(parameters)
  gradients[penaltyFunctionArguments$regularizedParameterLabels] <- tuningParameters$lambdas[penaltyFunctionArguments$regularizedParameterLabels]*
    parameters[penaltyFunctionArguments$regularizedParameterLabels]*
    (1/sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels])^2 + penaltyFunctionArguments$eps))
  return(gradients)
}

#' .smoothAdaptiveLASSOHessian
#' 
#' smoothed version of non-differentiable adaptive LASSO Hessian
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambdas (vector with one tuning parameter value for each parameter)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns Hessian matrix
#' @keywords internal
.smoothAdaptiveLASSOHessian <- function(parameters, 
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
#' .smoothElasticNetValue
#' 
#' smoothed version of non-differentiable elastic LASSO penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambda (tuning parameter value), alpha (0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns penalty function value
#' @keywords internal
.smoothElasticNetValue <- function(parameters, 
                             tuningParameters,
                             penaltyFunctionArguments){
  # using the implementation from https://github.com/psyphh/lslx/blob/ver.0.6.11/src/lslxOptimizer.cpp
  lassoPart <- .smoothLASSOValue(parameters = parameters,
                           penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                           "lambda" = tuningParameters$alpha*tuningParameters$lambda,
                                                           "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- .ridgeValue(parameters = parameters,
                     penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                     "lambda" = tuningParameters$lambda*(1-tuningParameters$alpha)
                     )
  )
  return(lassoPart + ridgePart)
}

#' .smoothElasticNetGradient
#' 
#' smoothed version of non-differentiable elastic LASSO gradient
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambda (tuning parameter value), alpha (0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns gradient values
#' @keywords internal
.smoothElasticNetGradient <- function(parameters, 
                                     tuningParameters,
                                     penaltyFunctionArguments){
  # using the implementation from https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
  lassoPart <- .smoothLASSOGradient(parameters = parameters,
                                   penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                                   "lambda" = tuningParameters$alpha*tuningParameters$lambda,
                                                                   "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- .ridgeGradient(parameters = parameters,
                             penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                             "lambda" = tuningParameters$lambda*(1-tuningParameters$alpha)
                             )
  )
  return(lassoPart + ridgePart)
}

#' .smoothElasticNetHessian
#' 
#' smoothed version of non-differentiable elastic LASSO Hessian
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with fields lambda (tuning parameter value), alpha (0<alpha<1. Controls the weight of ridge and lasso terms. alpha = 1 is lasso, alpha = 0 ridge)
#' @param penaltyFunctionArguments list with fields regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns Hessian matrix
#' @keywords internal
.smoothElasticNetHessian <- function(parameters, 
                                    tuningParameters,
                                    penaltyFunctionArguments){
  # using the implementation from https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
  lassoPart <- .smoothLASSOHessian(parameters = parameters,
                                  penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                                  "lambda" = tuningParameters$alpha*tuningParameters$lambda,
                                                                  "eps" = penaltyFunctionArguments$eps)
  )
  ridgePart <- .ridgeHessian(parameters = parameters,
                            penaltyFunctionArguments = list("regularizedParameterLabels" = penaltyFunctionArguments$regularizedParameterLabels,
                                                            "lambda" = tuningParameters$lambda*(1-tuningParameters$alpha)
                            )
  )
  return(lassoPart + ridgePart)
}

#### capped L1 ####
#' .smoothCappedL1Value
#' 
#' smoothed version of capped L1 penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns penalty function value
#' @keywords internal
.smoothCappedL1Value <- function(parameters, 
                             tuningParameters,
                             penaltyFunctionArguments
){
  smoothAbs <- sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels ])^2 + 
                      penaltyFunctionArguments$eps)
  penalty <- sum(tuningParameters$lambda*
                   ifelse(smoothAbs < tuningParameters$theta, smoothAbs, tuningParameters$theta)
                 )
  return(penalty)
}

#### SCAD ####
#' .smoothScadValue
#' 
#' smoothed version of scad penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns penalty function value
#' @keywords internal
.smoothScadValue <- function(parameters, 
                                tuningParameters,
                                penaltyFunctionArguments
){
  smoothAbs <- sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels ])^2 + 
                      penaltyFunctionArguments$eps)
  
  penalty <- 0
  
  for(p in 1:length(smoothAbs)){
    
    if(smoothAbs[p] <= tuningParameters$lambda){
      
      penalty <- penalty + tuningParameters$lambda * smoothAbs[p]
      
    }else if((tuningParameters$lambda < smoothAbs[p]) &&
             (smoothAbs[p] <= tuningParameters$lambda * tuningParameters$theta)
             ){
      
      penalty <- penalty + (- smoothAbs[p]^2 + 
                              2 * tuningParameters$theta * tuningParameters$lambda * 
                              smoothAbs[p] - 
                              tuningParameters$lambda^2) /(
                                2*(tuningParameters$theta-1)
                              ) 
      
    }else if(smoothAbs[p] > tuningParameters$lambda * tuningParameters$theta){
      
      penalty <- penalty + .5*(tuningParameters$theta + 1) * tuningParameters$lambda^2 
        
    }else{
      stop("impossible scad value")
    }
    
  }
  
  return(penalty)
}

#### LSP ####
#' .smoothLspValue
#' 
#' smoothed version of lsp penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns penalty function value
#' @keywords internal
.smoothLspValue <- function(parameters, 
                           tuningParameters,
                           penaltyFunctionArguments
){
  smoothAbs <- sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels ])^2 + 
                      penaltyFunctionArguments$eps)
  
  penalty <- tuningParameters$lambda * 
    log(1.0 + smoothAbs / tuningParameters$theta)
  
  
  return(sum(penalty))
}

#### MCP ####
#' .smoothMcpValue
#' 
#' smoothed version of mcp penalty
#' @param parameters vector with labeled parameter values
#' @param tuningParameters list with field lambda (tuning parameter value)
#' @param penaltyFunctionArguments list with field regularizedParameterLabels (labels of regularized parameters), and eps (controls the smooth approximation of non-differential penalty functions (e.g., lasso, adaptive lasso, or elastic net). Smaller values result in closer approximation, but may also cause larger issues in optimization.)
#' @returns penalty function value
#' @keywords internal
.smoothMcpValue <- function(parameters, 
                            tuningParameters,
                            penaltyFunctionArguments
){
  smoothAbs <- sqrt((parameters[penaltyFunctionArguments$regularizedParameterLabels ])^2 + 
                      penaltyFunctionArguments$eps)
  
  penalty <- 0
  
  for(p in 1:length(smoothAbs)){
    
    if(smoothAbs[p] <= (tuningParameters$theta * tuningParameters$lambda)){
      
      penalty <- penalty + tuningParameters$lambda * smoothAbs[p] -
        .5*(1/tuningParameters$theta)*smoothAbs[p]^2 
      
    }else if(smoothAbs[p] > (tuningParameters$theta * tuningParameters$lambda)
    ){
      
      penalty <- penalty + .5*(tuningParameters$lambda^2)*tuningParameters$theta
      
    }else{
      stop("impossible mcp value")
    }
    
  }
  
  return(penalty)
}