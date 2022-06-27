#' getMaxLambda_C
#' 
#' generates a the first lambda which sets all regularized parameters to zero
#' @param regularizedModel Model combining likelihood and lasso type penalty
#' @param SEM model of class Rcpp_SEMCpp
#' @param rawParameters labeled vector with starting values
#' @param weights weights given to each parameter in the penalty function
#' @param N sample size
getMaxLambda_C <- function(regularizedModel, 
                         SEM,
                         rawParameters,
                         weights,
                         N){
  
  lambda <- .Machine$double.xmax^(.05)
  result <- regularizedModel$optimize(
    rawParameters,
    SEM,
    lambda,
    1
  )
  
  sparseParameters <- result$rawParameters
  SEM <- linr:::setParameters(SEM = SEM, 
                              labels = names(sparseParameters), 
                              values = sparseParameters, 
                              raw = TRUE)
  SEM <- linr:::fit(SEM = SEM)
  gradients <- linr:::getGradients(SEM = SEM, 
                                   raw = TRUE)
  
  # define maxLambda as the maximal gradient of the regularized parameters
  maxLambda <- max(abs(gradients[weights != 0]) * 
                     weights[weights != 0]^(-1))
  # reset SEM
  SEM <- linr:::setParameters(SEM = SEM, 
                              labels = names(rawParameters), 
                              values = rawParameters, 
                              raw = TRUE)
  SEM <- linr:::fit(SEM = SEM)
  
  return((1/N)*(maxLambda+.1*maxLambda)) # adding some wiggle room as well
}