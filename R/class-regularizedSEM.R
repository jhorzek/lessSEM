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

setMethod("coef", "regularizedSEM", function (object, lambda = NULL, alpha = NULL) {
  if(is.null(lambda) && is.null(alpha)) return(object@parameters)
  if(is.null(lambda)) return(object@parameters[object@parameters$alpha == alpha,])
  if(is.null(alpha)) return(object@parameters[object@parameters$lambda == lambda,])
  
  pars <- object@parameters
  pars <- unlist(pars[pars$lambda == lambda & pars$alpha == alpha, object@parameterLabels])
  return(pars)
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
