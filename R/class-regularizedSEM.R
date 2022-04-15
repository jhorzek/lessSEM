setClass(Class = "regularizedSEM",
         representation = representation(
           parameters="data.frame",
           fits="data.frame", 
           internalOptimization="list", 
           inputArguments="list"
         )
)

setMethod("coef", "regularizedSEM", function (object, lambda = NULL, alpha = NULL) {
  if(is.null(lambda) && is.null(alpha)) return(object@parameters)
  if(any(is.null(lambda), is.null(alpha))) stop("Please specify both, lambda and alpha")
  if(any(!lambda %in% unique(object@parameters$lambda), !alpha %in% unique(object@parameters$alpha))
     ) stop("Could not find the requested alpha and lambda values in the model.")
  pars <- object@parameters
  pars <- pars[pars$lambda == lambda & pars$alpha == alpha, , drop = FALSE]
  parameters <- pars$value
  names(parameters) <- pars$label
  return(parameters)
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
    parameters <- parameters[parameters$alpha == alpha & parameters$regularized, , drop = FALSE]
    
    ggplot2::ggplot(data = parameters,
                                     mapping = ggplot2::aes(x = lambda, y = value, group = label)) +
      ggplot2::geom_line(colour = "#008080")+
      ggplot2::ggtitle("Regularized Parameters")
    
  }else{
    parameters <- parameters[parameters$alpha == alpha,]
    
    ggplot2::ggplot(data = parameters,
                                     mapping = ggplot2::aes(x = lambda, y = value, group = label, color = regularized)) +
      ggplot2::geom_line()+ 
      ggplot2::scale_color_manual(values=c("FALSE"="gray","TRUE"="#008080")) +
      ggplot2::ggtitle("Parameter Estimates")
  }
})
