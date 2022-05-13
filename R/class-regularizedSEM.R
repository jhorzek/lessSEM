setClass(Class = "regularizedSEM",
         representation = representation(
           parameters="data.frame",
           fits="data.frame", 
           parameterLabels = "character",
           regularized = "character",
           internalOptimization="list", 
           inputArguments="list"
         )
)

#' show
#' @param object object of class regularizedSEM
#' @export
setMethod("show", "regularizedSEM", function (object) {
  #modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class regularizedSEM with ",object@inputArguments$penalty, " penalty ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat(paste0("- Use coef(object) to get the parameter estimates of the model. With coef(object, lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Use plot(object) to plot the parameter estimates of the model.\n\n"))
  cat(paste0("- Use aCV4regularizedSEM(object, k = k) to compute an approximate k-fold cross-valdidation.\n\n"))
  cat(paste0("- Information criteria can be compute with AIC(object) or BIC(object).\n\n"))
  cat("################################################\n")
})

#' summary
#' @param object object of class regularizedSEM
#' @export
setMethod("summary", "regularizedSEM", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class regularizedSEM with ",object@inputArguments$penalty, " penalty ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat(paste0("- Use coef(", modelName, 
             ") to get the parameter estimates of the model. With coef(", 
             modelName, "lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Use plot(", modelName, 
             ") to plot the parameter estimates of the model.\n\n"))
  cat(paste0("- Use aCV4regularizedSEM(", modelName, 
             ", k = k) to compute an approximate k-fold cross-valdidation.\n\n"))
  cat(paste0("- Information criteria can be compute with AIC(", modelName, 
             ") or BIC(", modelName, 
             ").\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of a regularizedSEM
#' 
#' @param object object of class regularizedSEM
#' @param criterion can be one of: "AIC", "BIC". If set to NULL, all parameters will be returned
#' @param lambda numeric: value of lambda for which parameters should be returned. Must be present in object
#' @param alpha numeric: value of alpha for which parameters should be returned. Must be present in object
#' @export
setMethod("coef", "regularizedSEM", function (object, criterion = NULL, lambda = NULL, alpha = NULL) {
  if(is.null(lambda) && is.null(alpha) && is.null(criterion)) return(object@parameters)
  if(is.null(lambda) && is.null(criterion)) return(object@parameters[object@parameters$alpha == alpha,])
  if(is.null(alpha) && is.null(criterion)) return(object@parameters[object@parameters$lambda == lambda,])
  if(!is.null(criterion) && criterion %in% c("AIC", "BIC")){
    if(length(unique(object@parameters$alpha)) != 1 || unique(object@parameters$alpha) != 1) stop("Selection by criterion currently only supported for lasso type penalties.")
    if(criterion == "AIC"){
      AICs <- AIC(object)
      bestAIC <- which(AICs$AIC == min(AICs$AIC))[1]
      return(object@parameters[bestAIC,]) 
    }
    
    if(criterion == "BIC"){
      BICs <- BIC(object)
      bestBIC <- which(BICs$BIC == min(BICs$BIC))[1]
      return(object@parameters[bestBIC,])
    }
    
  }
  
  pars <- object@parameters
  pars <- unlist(pars[pars$lambda == lambda & pars$alpha == alpha, object@parameterLabels])
  return(pars)
})

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class regularizedSEM
#' @export
setMethod("AIC", "regularizedSEM", function (object) {
  fits <- object@fits
  if(all(fits$alpha == 1)){
    fits$AIC <- fits$m2LL + 2*fits$nonZeroParameters
  }else{
    warning("AIC for non-lasso models is experimental and should not be trusted. Use approximate cross-validation instead")
    fits$AIC <- rep(NA, length(fits$m2LL))
    tuningParameters <- data.frame(lambda = fits$lambda, alpha = fits$alpha)
    
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = object@inputArguments$lavaanModel, fit = FALSE)
    
    for(i in 1:nrow(tuningParameters)){
      SEM <- aCV4SEM:::setParameters(SEM, 
                                     labels = object@parameterLabels, 
                                     values = coef(object, 
                                                   lambda = tuningParameters$lambda[i], 
                                                   alpha = tuningParameters$alpha[i]), 
                                     raw = FALSE)
      SEM <- aCV4SEM:::fit(SEM)
      scores <- aCV4SEM:::getScores(SEM, raw = FALSE)
      hessian <- aCV4SEM:::getHessian(SEM, raw = FALSE)
      hessianInv <- solve(hessian)
      twoTimesNpar <- sum(apply(scores, 1, function(x) t(x)%*%hessianInv%*%(x)))
      fits$AIC[i] <- fits$m2LL[i] + twoTimesNpar
    }
  }
  
  return(fits)
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class regularizedSEM
#' @export
setMethod("BIC", "regularizedSEM", function (object) {
  N <- nrow(lavaan::lavInspect(object@inputArguments$lavaanModel, "data"))
  fits <- object@fits
  
  if(all(fits$alpha == 1)){
    fits$BIC <- fits$m2LL + log(N)*fits$nonZeroParameters
  }else{
    warning("BIC for non-lasso models is experimental and should not be trusted. Use approximate cross-validation instead")
    fits$BIC <- rep(NA, length(fits$m2LL))
    tuningParameters <- data.frame(lambda = fits$lambda, alpha = fits$alpha)
    
    SEM <- aCV4SEM:::SEMFromLavaan(lavaanModel = object@inputArguments$lavaanModel, fit = FALSE)
    
    for(i in 1:nrow(tuningParameters)){
      SEM <- aCV4SEM:::setParameters(SEM, 
                                     labels = object@parameterLabels, 
                                     values = coef(object, 
                                                   lambda = tuningParameters$lambda[i], 
                                                   alpha = tuningParameters$alpha[i]), 
                                     raw = FALSE)
      SEM <- aCV4SEM:::fit(SEM)
      scores <- aCV4SEM:::getScores(SEM, raw = FALSE)
      hessian <- aCV4SEM:::getHessian(SEM, raw = FALSE)
      hessianInv <- solve(hessian)
      npar <- .5*sum(apply(scores, 1, function(x) t(x)%*%hessianInv%*%(x))) # scores and hessian are for -2 log-Likelihood
      fits$BIC[i] <- fits$m2LL[i] + log(N)*npar
    }
  }

  return(fits)
})

#' plot
#' 
#' plots the regularized and unregularized parameters for all levels of lambda
#' 
#' @param x object of class regularizedSEM
#' @param alpha numeric: value of alpha for which parameters should be returned. Must be present in object. Required if elastic net was used
#' @param regularizedOnly boolean: should only regularized parameters be plotted?``
#' @export
setMethod("plot", "regularizedSEM", function (x, alpha = NULL, regularizedOnly = TRUE) {
  parameters <- x@parameters
  if(is.null(alpha)) {
    alpha <- unique(parameters$alpha)[1]
    
    if(length(unique(parameters$alpha)) != 1) warning(paste0("Models for different values of alpha were fitted, but not alpha was provided to the plotting function. Plotting for alpha=", 
                                                             alpha, ". You can specify specific alpha values with, for instance, plot(aCV4Regsem, alpha = .02)"))
  }
  if(!alpha %in% parameters$alpha){
    stop(paste0("alpha=", alpha, " is was not tested in the model. Only found the following alpha values: "),
         paste0(unique(parameters$alpha), collapse= ", "), "."
    )
  }
  
  if(regularizedOnly){
    parameters <- cbind(
      lambda = x@fits$lambda,
      alpha = x@fits$alpha,
      coef(x, alpha = alpha)[,x@regularized, drop = FALSE]
    )
    parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@regularized)
    
    ggplot2::ggplot(data = parametersLong,
                    mapping = ggplot2::aes(x = lambda, y = value, group = name)) +
      ggplot2::geom_line(colour = "#008080")+
      ggplot2::ggtitle("Regularized Parameters")
    
  }else{
    parameters <- cbind(
      lambda = x@fits$lambda,
      alpha = x@fits$alpha,
      coef(x, alpha = alpha)
    )
    parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@parameterLabels)
    parametersLong$regularized <- parametersLong$name %in% x@regularized
    
    ggplot2::ggplot(data = parametersLong,
                    mapping = ggplot2::aes(x = lambda, 
                                           y = value, 
                                           group = name, 
                                           color = regularized)) +
      ggplot2::geom_line()+ 
      ggplot2::scale_color_manual(values=c("FALSE"="gray","TRUE"="#008080")) +
      ggplot2::ggtitle("Parameter Estimates")
  }
})
