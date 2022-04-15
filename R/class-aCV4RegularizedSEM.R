setClass(Class = "aCV4RegularizedSEM",
         representation = representation(
           parameters="data.frame",
           cvfits = "data.frame",
           cvfitsDetails="data.frame", 
           subsets = "list"
         )
)

setMethod("show", "aCV4RegularizedSEM", function (object) {
  bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
  cat(paste0("Best parameter estimates observed for: lambda = ", object@cvfits$lambda[bestFit], " & alpha = ", object@cvfits$alpha[bestFit], ":\n"))
  
  parameters <- object@parameters$value[object@parameters$id == object@cvfits$id[bestFit]]
  names(parameters) <- object@parameters$label[object@parameters$id == object@cvfits$id[bestFit]]
  print(parameters)
})

setMethod("coef", "aCV4RegularizedSEM", function (object, lambda = NULL, alpha = NULL) {
  
  if(is.null(lambda) && is.null(alpha)) {
    # return parameters of model with best cv fit
    bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
    parameters <- object@parameters$value[object@parameters$id == object@cvfits$id[bestFit]]
    names(parameters) <- object@parameters$label[object@parameters$id == object@cvfits$id[bestFit]]
    return(parameters)
  }
  # otherwise: return the parameters for the requested lambda and alpha values:
  if(any(is.null(lambda), is.null(alpha))) stop("Please specify both, lambda and alpha")
  if(any(!lambda %in% unique(object@parameters$lambda), !alpha %in% unique(object@parameters$alpha))
  ) stop("Could not find the requested alpha and lambda values in the model.")
  pars <- object@parameters
  pars <- pars[pars$lambda == lambda & pars$alpha == alpha, , drop = FALSE]
  parameters <- pars$value
  names(parameters) <- pars$label
  return(parameters)
  
})


setMethod("plot", "aCV4RegularizedSEM", function (x, alpha = NULL, regularizedOnly = TRUE) {
  parameters <- x@parameters
  fits <- x@cvfits
  if(is.null(alpha)) {
    alpha <- unique(parameters$alpha)[1]
    
    if(length(unique(parameters$alpha)) != 1) warning(paste0("Models for different values of alpha were fitted, but not alpha was provided to the plotting function. Plotting for alpha=", 
                                                             alpha, ". You can specify specific alpha values with, for instance, plot(aCV4RegularizedSEM, alpha = .02)"))
  }
  if(!alpha %in% parameters$alpha){
    stop(paste0("alpha=", alpha, " is was not tested in the model. Only found the following alpha values: "),
         paste0(unique(parameters$alpha), collapse= ", "), "."
    )
  }
  
  if(regularizedOnly){
    parameters <- parameters[parameters$alpha == alpha & parameters$regularized, , drop = FALSE]
    fits <- fits[fits$alpha == alpha,, drop = FALSE]
    
    parameterPlot <- ggplot2::ggplot(data = parameters,
                                     mapping = ggplot2::aes(x = lambda, y = value, group = label)) +
      ggplot2::geom_line(colour = "#008080")+
      ggplot2::ggtitle("Regularized Parameters")
    fitPlot <- ggplot2::ggplot(data = fits,
                               mapping = ggplot2::aes(x = lambda, y = cvfit), title = "cv-fit") +
      ggplot2::geom_line()

  }else{
    parameters <- parameters[parameters$alpha == alpha,]
    fits <- fits[fits$alpha == alpha,]
    
    parameterPlot <- ggplot2::ggplot(data = parameters,
                                     mapping = ggplot2::aes(x = lambda, y = value, group = label, color = regularized)) +
      ggplot2::geom_line()+ 
      ggplot2::scale_color_manual(values=c("FALSE"="gray","TRUE"="#008080")) +
      ggplot2::ggtitle("Parameter Estimates")
    fitPlot <- ggplot2::ggplot(data = fits,
                               mapping = ggplot2::aes(x = lambda, y = cvfit), title = "cv-fit") +
      ggplot2::geom_line()

  }

  parameterPlot / fitPlot
})