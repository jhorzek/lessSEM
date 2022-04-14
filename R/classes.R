# https://gallery.rcpp.org/articles/custom-printer-exposed-modules/
setClass("Rcpp_SEMCpp")
setMethod("show", "Rcpp_SEMCpp", function (object) {
  cat("Internal C++ model representation of aCV4SEM\n")
  cat("Parameters:\n")
  print(aCV4SEM::getParameters(object))
  cat("\n")
  cat(paste0("-2 log-Likelihood: ", object$m2LL))
})

setMethod("logLik", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  N <- nrow(object$rawData)
  numberOfParameters <- length(getParameters(object))
  
  ll <- new("logLikelihood",
            logLik = -.5*object$m2LL,
            nParameters = numberOfParameters,
            N = N)
  return(ll)
})

setClass("logLikelihood",
         representation = representation(
           logLik="numeric",
           nParameters="integer", 
           N="integer"
         ))
setMethod("show", "logLikelihood", function(object){
  cat(paste0(object@logLik), "(nPar = ", object@nParameters, ")")
})

setMethod("coef", "Rcpp_SEMCpp", function (object) {
  return(getParameters(object, raw = FALSE))
})

setMethod("AIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  AICis <- -2*ll@logLik + 2*ll@nParameters
  return(AICis)
})

setMethod("BIC", "Rcpp_SEMCpp", function (object) {
  if(!object$wasFit){
    object$fit()
  }
  ll <- logLik(object)
  BICis <- -2*ll@logLik + log(ll@N)*ll@nParameters
  return(BICis)
})

setClass(Class = "regularizedSEM",
         representation = representation(
           parameters="data.frame",
           fits="data.frame", 
           internalOptimization="list", 
           inputArguments="list"
         )
)

setMethod("coef", "regularizedSEM", function (object) {
  return(object@parameters)
})

setMethod("AIC", "regularizedSEM", function (object) {
  fits <- object@fits
  fits$AIC <- fits$m2LL + 2*fits$nonZeroParameters
  
  return(fits)
})

setMethod("BIC", "regularizedSEM", function (object) {
  N <- nrow(lavaan::lavInspect(object@inputArguments$lavaanModel, "data"))
  fits <- object@fits
  fits$BIC <- fits$m2LL + log(N)*fits$nonZeroParameters
  
  return(fits)
})

setMethod("plot", "regularizedSEM", function (x, alpha = NULL) {
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
  
  parameters <- parameters[parameters$alpha == alpha,]
  ggplot2::ggplot(data = parameters, mapping = ggplot2::aes(x = lambda, y = value, group = label, color = regularized))+
    ggplot2::geom_line()
})

setClass(Class = "aCV4Regsem",
         representation = representation(
           parameters="data.frame",
           cvfits = "data.frame",
           cvfitsDetails="data.frame", 
           subsets = "list"
         )
)

setMethod("show", "aCV4Regsem", function (object) {
  bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
  cat(paste0("Best parameter estimates observed for: lambda = ", object@cvfits$lambda[bestFit], " & alpha = ", object@cvfits$alpha[bestFit], ":\n"))
  
  parameters <- object@parameters$value[object@parameters$id == object@cvfits$id[bestFit]]
  names(parameters) <- object@parameters$label[object@parameters$id == object@cvfits$id[bestFit]]
  print(parameters)
})

setMethod("coef", "aCV4Regsem", function (object) {
  bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
  parameters <- object@parameters$value[object@parameters$id == object@cvfits$id[bestFit]]
  names(parameters) <- object@parameters$label[object@parameters$id == object@cvfits$id[bestFit]]
  return(parameters)
})


setMethod("plot", "aCV4Regsem", function (x, alpha = NULL) {
  parameters <- x@parameters
  fits <- x@cvfits
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
  
  parameters <- parameters[parameters$alpha == alpha,]
  fits <- fits[fits$alpha == alpha,]
  
  parameterPlot <- ggplot2::ggplot(data = parameters,
                                   mapping = ggplot2::aes(x = lambda, y = value, group = label, color = regularized))+
    ggplot2::geom_line()
  fitPlot <- ggplot2::ggplot(data = fits,
                             mapping = ggplot2::aes(x = lambda, y = cvfit), title = "cv-fit")+
    ggplot2::geom_line()
  parameterPlot / fitPlot
})