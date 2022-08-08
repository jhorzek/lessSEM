#' GIC
#' 
#' THIS FUNCTION IS UNDER DEVELOPMENT AND SHOULD _NOT_ BE USED.
#' computes the generalized information criterion as 
#' \eqn{-2\log \mathcal L + \text{scaler}*\text{df}}
#' where scaler is a numeric value with which the 
#' degrees of freedom (df) are multiplied.
#' 
#' See Zhang et al. (2010) for more details.
#' 
#' * Zhang, Y., Li, R., & Tsai, C.-L. (2010). Regularization Parameter 
#' Selections via Generalized Information Criterion. 
#' Journal of the American Statistical Association, 105(489), 
#' 312–323. https://doi.org/10.1198/jasa.2009.tm08013
#' @param regularizedSEM model of class regularizedSEM
#' @param scaler numeric value to scale the degrees of freedom
#' @return vector with GIC values.
GIC <- function(regularizedSEM, scaler = 2){
  warning("GIC IS EXPERIMENTAL AND SHOULD NOT BE USED!")
  penalty <- regularizedSEM@penalty
  
  if(!penalty %in% c("scad", "mcp", "adaptiveLasso")) message("GIC has been proposed for consistent estimators.")
  
  # check if regularizedSEM uses approximated penalty function
  if(!is.null(regularizedSEM@inputArguments$epsilon)){
    epsilon <- regularizedSEM@inputArguments$epsilon
  }else{
    epsilon <- 0
  }
  
  # compute Hessian for each value of the tuning parameters
  parameters <- regularizedSEM@parameters[,regularizedSEM@parameterLabels]
  
  gic <- rep(NA, nrow(parameters))
  
  # we need a model to compute the Hessians
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = regularizedSEM@inputArguments$lavaanModel, 
                                whichPars = "start", 
                                transformVariances = TRUE, 
                                fit = FALSE,
                                addMeans = regularizedSEM@inputArguments$modifyModel$addMeans,
                                activeSet = regularizedSEM@inputArguments$modifyModel$activeSet, 
                                dataSet = regularizedSEM@inputArguments$modifyModel$dataSet
  )
  
  N <- nrow(SEM$rawData)
  
  pbar <- utils::txtProgressBar(min = 0, max = nrow(parameters), style = 3)
  
  for(p in 1:nrow(parameters)){
    
    utils::setTxtProgressBar(pb = pbar, value = p)
    
    SEM <- lessSEM:::.setParameters(SEM = SEM, 
                         labels = regularizedSEM@parameterLabels, 
                         values = unlist(parameters[p,]), 
                         raw = FALSE)
    
    m2LL <- lessSEM:::.fit(SEM = SEM)$m2LL
    
    # We need the Hessian of the -2log-likelihood and the penalty function
    m2LLHessian <- lessSEM:::.getHessian(SEM, raw = FALSE)
    
    if(penalty == "lasso"){
      
      penaltyHessian <- N*smoothLASSOHessian(parameters = unlist(parameters[p,]), 
                                             tuningParameters = list(
                                               lambda = regularizedSEM@fits$lambda[p]
                                             ), 
                                             penaltyFunctionArguments = list(
                                               regularizedParameterLabels = regularizedSEM@regularized,
                                               eps = epsilon
                                             )
      )
      
    }
    
    # remove all rows and columns with NAs; these are the ones, where parameters 
    # have been zeroed
    rowWithNA <- apply(penaltyHessian,1,anyNA)
    colWithNA <- apply(penaltyHessian,2,anyNA)
    
    penaltyHessian <- penaltyHessian[!rowWithNA,]
    penaltyHessian <- penaltyHessian[,!colWithNA]
    
    m2LLHessian <- m2LLHessian[!rowWithNA,]
    m2LLHessian <- m2LLHessian[,!colWithNA]
    
    dfs <- try(sum(diag(solve(m2LLHessian - penaltyHessian)%*%m2LLHessian)))
    if(is(dfs, "try-error")) next
    
    gic[p] <- m2LL + scaler*dfs
  }
  
  return(gic)
}
