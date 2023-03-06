#' GIC
#' 
#' THIS FUNCTION IS UNDER DEVELOPMENT AND SHOULD _NOT_ BE USED.
#' computes the generalized information criterion as 
#' 
#' -2-log-Likelihood + k*df
#' 
#' where the k is a numeric value with which the 
#' degrees of freedom (df) are multiplied. To get the equivalent
#' of the AIC, use k = 2. For an equivalent to the BIC, set k = log(N)
#' 
#' See Fan & Li (2001), p. 1355 and Zhang et al. (2010), p. 314 for more details.
#' 
#' * Zhang, Y., Li, R., & Tsai, C.-L. (2010). Regularization Parameter 
#' Selections via Generalized Information Criterion. 
#' Journal of the American Statistical Association, 105(489), 
#' 312–323. https://doi.org/10.1198/jasa.2009.tm08013
#' @param regularizedSEM model of class regularizedSEM
#' @param k numeric value to scale the degrees of freedom
#' @return vector with GIC values.
GIC <- function(regularizedSEM, k = 2){
  warning("GIC IS EXPERIMENTAL AND SHOULD NOT BE USED!")
  penalty <- regularizedSEM@penalty
  
  if(!penalty %in% c("scad", "mcp", "adaptiveLasso")) .printNote("GIC has been proposed for consistent estimators.")
  
  # check if regularizedSEM uses approximated penalty function
  if(!is.null(regularizedSEM@inputArguments$epsilon)){
    epsilon <- regularizedSEM@inputArguments$epsilon
  }else{
    epsilon <- 0
  }
  
  # compute Hessian for each value of the tuning parameters
  parameters <- regularizedSEM@parameters[,regularizedSEM@parameterLabels]
  
  gic <- rep(NA, nrow(parameters))
  dfs <- gic
  
  # we need a model to compute the Hessians
  SEM <- .SEMFromLavaan(lavaanModel = regularizedSEM@inputArguments$lavaanModel, 
                        whichPars = "start",
                        fit = FALSE,
                        addMeans = regularizedSEM@inputArguments$modifyModel$addMeans,
                        activeSet = regularizedSEM@inputArguments$modifyModel$activeSet, 
                        dataSet = regularizedSEM@inputArguments$modifyModel$dataSet,
                        transformations = regularizedSEM@inputArguments$modifyModel$transformations,
                        transformationList = regularizedSEM@inputArguments$modifyModel$transformationList
  )
  
  N <- nrow(SEM$rawData)
  
  pbar <- utils::txtProgressBar(min = 0, max = nrow(parameters), style = 3)
  
  for(p in 1:nrow(parameters)){
    
    utils::setTxtProgressBar(pb = pbar, value = p)
    
    SEM <- .setParameters(SEM = SEM, 
                          labels = regularizedSEM@parameterLabels, 
                          values = unlist(parameters[p,]), 
                          raw = FALSE)
    
    m2LL <- .fit(SEM = SEM)$m2LL
    
    # We need the Hessian of the -2log-likelihood and the penalty function
    m2LLHessian <- .getHessian(SEM, raw = TRUE)
    
    if(penalty == "lasso"){
      
      penaltyHessian <- N*.smoothLASSOHessian(parameters = unlist(parameters[p,]), 
                                             tuningParameters = list(
                                               lambda = regularizedSEM@fits$lambda[p]
                                             ), 
                                             penaltyFunctionArguments = list(
                                               regularizedParameterLabels = regularizedSEM@regularized,
                                               eps = epsilon
                                             )
      )
      
    }else{
      stop("Currently only implemented for lasso penalty")
    }
    
    # remove all rows and columns with NAs; these are the ones, where parameters 
    # have been zeroed
    rowWithNA <- apply(penaltyHessian,1,anyNA)
    colWithNA <- apply(penaltyHessian,2,anyNA)
    
    penaltyHessian <- penaltyHessian[!rowWithNA,]
    penaltyHessian <- penaltyHessian[,!colWithNA]
    
    m2LLHessian <- m2LLHessian[!rowWithNA,]
    m2LLHessian <- m2LLHessian[,!colWithNA]
    
    # Note: the fomula in Zhang et al. (2010), p. 314 is based on the log-
    # likelihood, note the -2-log-likelihood -> multiplication with -.5
    df <- try(sum(diag(solve(-.5*m2LLHessian - .5*penaltyHessian)%*%(-.5*m2LLHessian))))
    if(is(df, "try-error")) next
    dfs[p] <- df
    gic[p] <- m2LL + k*df
  }
  
  return(
    list("gic" = gic,
         "df" = dfs)
  )
}
