setClass(Class = "regularizedSEMWithCustomPenalty",
         representation = representation(
           parameters="data.frame",
           fits="data.frame", 
           parameterLabels = "character",
           inputArguments="list"
         )
)

#' summary
#' @param object object of class regularizedSEMWithCustomPenalty
#' @export
setMethod("summary", "regularizedSEMWithCustomPenalty", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class regularizedSEMWithCustomPenalty ####\n\n"))
  cat(paste0("- Use coef(", modelName, 
             ") to get the parameter estimates of the model. With coef(", 
             modelName, "lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Use plot(", modelName, 
             ") to plot the parameter estimates of the model.\n\n"))
  cat(paste0("- Information criteria can be compute with AIC(", modelName, 
             ") or BIC(", modelName, 
             ").\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of a regularizedSEMWithCustomPenalty
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @export
setMethod("coef", "regularizedSEMWithCustomPenalty", function (object) {
  pars <- object@parameters
  return(pars)
})

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @export
setMethod("AIC", "regularizedSEMWithCustomPenalty", function (object) {
  stop("AIC is experimental and should not be trusted. Use cross-validation instead, Not correctly implemented yet.")
  fits <- object@fits
  fits$AIC <- rep(NA, length(fits$m2LL))
  tuningParameters <- object@inputArguments$tuningParameters
  penaltyFunctionArguments <- object@inputArguments$penaltyFunctionArguments
  
  SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = object@inputArguments$lavaanModel, fit = FALSE)
  N <- nrow(SEM$rawData)
  
  for(i in 1:nrow(tuningParameters)){
    parameters <- unlist(object@parameters[i,object@parameterLabels])
    SEM <- aCV4SEM:::setParameters(SEM, 
                                   labels = object@parameterLabels, 
                                   values = parameters,
                                   raw = FALSE)
    SEM <- aCV4SEM:::fit(SEM)
    scores <- aCV4SEM:::getScores(SEM, raw = FALSE)
    hessian <- aCV4SEM:::getHessian(SEM, raw = FALSE)
    
    # penalty scores
    penaltyScores <- object@inputArguments$individualPenaltyFunctionGradient(parameters = parameters, 
                                                                             tuningParameters = tuningParameters[i,,drop=FALSE], 
                                                                             penaltyFunctionArguments = penaltyFunctionArguments)
    penaltyHessian <- N*object@inputArguments$individualPenaltyFunctionHessian(parameters = parameters, 
                                                                               tuningParameters = tuningParameters[i,,drop=FALSE], 
                                                                               penaltyFunctionArguments = penaltyFunctionArguments)
    
    scores <- scores + matrix(rep(penaltyScores, N), 
                              nrow = N, 
                              ncol = length(penaltyScores), 
                              byrow = TRUE,
                              dimnames = list(NULL, names(penaltyScores)))[,colnames(scores)]
    hessian <- hessian + penaltyHessian[rownames(hessian), colnames(hessian)]
    
    hessianInv <- solve(hessian)
    twoTimesNpar <- sum(apply(scores, 1, function(x) t(x)%*%hessianInv%*%(x)))
    fits$AIC[i] <- fits$m2LL[i] + twoTimesNpar
  }
  
  fits <- cbind(tuningParameters, fits)
  return(fits)
  
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class regularizedSEMWithCustomPenalty
#' @export
setMethod("BIC", "regularizedSEMWithCustomPenalty", function (object) {
  stop("BIC is experimental and should not be trusted. Use cross-validation instead. Not correctly implemented yet.")
  fits <- object@fits
  fits$BIC <- rep(NA, length(fits$m2LL))
  tuningParameters <- object@inputArguments$tuningParameters
  penaltyFunctionArguments <- object@inputArguments$penaltyFunctionArguments
  
  SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = object@inputArguments$lavaanModel, fit = FALSE)
  N <- nrow(SEM$rawData)
  
  for(i in 1:nrow(tuningParameters)){
    parameters <- unlist(object@parameters[i,object@parameterLabels])
    SEM <- aCV4SEM:::setParameters(SEM, 
                                   labels = object@parameterLabels, 
                                   values = parameters,
                                   raw = FALSE)
    SEM <- aCV4SEM:::fit(SEM)
    scores <- aCV4SEM:::getScores(SEM, raw = FALSE)
    hessian <- aCV4SEM:::getHessian(SEM, raw = FALSE)
    
    # penalty scores
    penaltyScores <- object@inputArguments$individualPenaltyFunctionGradient(parameters = parameters, 
                                                                             tuningParameters = tuningParameters[i,,drop=FALSE], 
                                                                             penaltyFunctionArguments = penaltyFunctionArguments)
    penaltyHessian <- N*object@inputArguments$individualPenaltyFunctionHessian(parameters = parameters, 
                                                                               tuningParameters = tuningParameters[i,,drop=FALSE], 
                                                                               penaltyFunctionArguments = penaltyFunctionArguments)
    
    scores <- scores + matrix(rep(penaltyScores, N), 
                              nrow = N, 
                              ncol = length(penaltyScores), 
                              byrow = TRUE,
                              dimnames = list(NULL, names(penaltyScores)))[,colnames(scores)]
    hessian <- hessian + penaltyHessian[rownames(hessian), colnames(hessian)]
    
    hessianInv <- solve(hessian)
    twoTimesNpar <- sum(apply(scores, 1, function(x) t(x)%*%hessianInv%*%(x)))
    fits$BIC[i] <- fits$m2LL[i] + log(N)*.5*twoTimesNpar
  }
  
  fits <- cbind(tuningParameters, fits)
  return(fits)
})
