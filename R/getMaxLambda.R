#' .getMaxLambda_C
#' 
#' generates a the first lambda which sets all regularized parameters to zero
#' @param regularizedModel Model combining likelihood and lasso type penalty
#' @param SEM model of class Rcpp_SEMCpp
#' @param rawParameters labeled vector with starting values
#' @param weights weights given to each parameter in the penalty function
#' @param N sample size
#' @param approx When set to TRUE, .Machine$double.xmax^(.01) is used instead of .Machine$double.xmax^(.05)
#' @returns first lambda value which sets all regularized parameters to zero (plus some tolerance)
#' @keywords internal
.getMaxLambda_C <- function(regularizedModel, 
                         SEM,
                         rawParameters,
                         weights,
                         N,
                         approx  = FALSE){
  
  lambda <- ifelse(approx,
                   .Machine$double.xmax^(.01),
                   .Machine$double.xmax^(.05)
  )
  result <- regularizedModel$optimize(
    rawParameters,
    SEM,
    lambda,
    1 # alpha = 1
  )
  
  sparseParameters <- result$rawParameters
  SEM <- .setParameters(SEM = SEM, 
                              labels = names(sparseParameters), 
                              values = sparseParameters, 
                              raw = TRUE)
  SEM <- .fit(SEM = SEM)
  gradients <- .getGradients(SEM = SEM, 
                                   raw = TRUE)
  
  # define maxLambda as the maximal gradient of the regularized parameters
  maxLambda <- max(abs(gradients[weights != 0]) * 
                     weights[weights != 0]^(-1))
  # reset SEM
  SEM <- .setParameters(SEM = SEM, 
                              labels = names(rawParameters), 
                              values = rawParameters, 
                              raw = TRUE)
  SEM <- .fit(SEM = SEM)
  
  return((1/N)*(maxLambda+.1*maxLambda)) # adding some wiggle room as well
}

#' .gpGetMaxLambda
#' 
#' generates a the first lambda which sets all regularized parameters to zero
#' @param regularizedModel Model combining likelihood and lasso type penalty
#' @param par labeled vector with starting values
#' @param fitFunction R fit function 
#' @param gradientFunction R gradient functions
#' @param userSuppliedArguments list with arguments for fitFunction and gradientFunction
#' @param weights weights given to each parameter in the penalty function
#' @returns first lambda value which sets all regularized parameters to zero (plus some tolerance)
#' @keywords internal
.gpGetMaxLambda <- function(regularizedModel,
                           par,
                           fitFunction,
                           gradientFunction,
                           userSuppliedArguments,
                           weights){
  
  lambda <- .Machine$double.xmax^(.05)
  result <- regularizedModel$optimize(
    par,
    fitFunction,
    gradientFunction,
    userSuppliedArguments,
    lambda,
    1
  )
  
  sparseParameters <- result$rawParameters
  gradients <- gradientFunction(sparseParameters,
                                names(sparseParameters),
                                userSuppliedArguments)
  
  # define maxLambda as the maximal gradient of the regularized parameters
  maxLambda <- max(abs(gradients[weights[names(par)] != 0]) * 
                     weights[weights[names(par)] != 0]^(-1))
  
  return(maxLambda+.1*maxLambda) # adding some wiggle room as well
}