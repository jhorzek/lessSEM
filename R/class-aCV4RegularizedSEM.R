setClass(Class = "aCV4RegularizedSEM",
         representation = representation(
           parameters="data.frame",
           cvfits = "data.frame",
           parameterLabels = "character",
           regularized = "character",
           cvfitsDetails="data.frame", 
           subsets = "list"
         )
)

#' show
#' 
#' @export
setMethod("show", "aCV4RegularizedSEM", function (object) {
  bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
  cat(paste0("Best parameter estimates observed for: lambda = ", object@cvfits$lambda[bestFit], " & alpha = ", object@cvfits$alpha[bestFit], ":\n"))
  
  parameters <- unlist(object@parameters[bestFit,object@parameterLabels])
  print(parameters)
})

#' summary
#' 
#' @export
setMethod("summary", "aCV4RegularizedSEM", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Approximate Cross Validation Results ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat("Final parameter estimates based on best out-of-sample fit:\n")
  print(coef(object))
  cat("\n\n")
  cat(paste0("- Use coef(", modelName, 
             ") to get the final parameter estimates of the model. With coef(", 
             modelName, "lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Use plot(", modelName, 
             ") to plot the parameter estimates and the cross-validation fit of the model.\n\n"))
  cat(paste0("- Cross-validation fits can be accessed with ", modelName, 
             "@cvfits.\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of an aCV4RegularizedSEM
#' 
#' @param object object of class regularizedSEM
#' @param lambda numeric: value of lambda for which parameters should be returned. Must be present in object
#' @param alpha numeric: value of alpha for which parameters should be returned. Must be present in object
#' 
#' @export
setMethod("coef", "aCV4RegularizedSEM", function (object, lambda = NULL, alpha = NULL) {
  
  if(is.null(lambda) && is.null(alpha)) {
    # return parameters of model with best cv fit
    bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
    parameters <- unlist(object@parameters[bestFit,object@parameterLabels])
    return(parameters)
  }
  
  if(is.null(lambda)) return(object@parameters[object@parameters$alpha == alpha,object@parameterLabels])
  if(is.null(alpha)) return(object@parameters[object@parameters$lambda == lambda,object@parameterLabels])
  
  pars <- object@parameters
  pars <- unlist(pars[pars$lambda == lambda & pars$alpha == alpha, object@parameterLabels])
  return(pars)
})

#' plot
#' 
#' plots the regularized and unregularized parameters as well as the cross-validation fits for all levels of lambda
#' 
#' @param x object of class regularizedSEM
#' @param alpha numeric: value of alpha for which parameters should be returned. Must be present in object. Required if elastic net was used
#' @param regularizedOnly boolean: should only regularized parameters be plotted?``
#' @export
setMethod("plot", "aCV4RegularizedSEM", function (x, alpha = NULL, regularizedOnly = TRUE) {
  
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
  
  if(regularizedOnly){
    fits <- fits[fits$alpha == alpha,, drop = FALSE]
    parameters <- cbind(
      lambda = parameters$lambda,
      alpha = parameters$alpha,
      coef(x, alpha = alpha)[,x@regularized, drop = FALSE]
    )
    parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@regularized)
    
    parameterPlot <- ggplot2::ggplot(data = parametersLong,
                    mapping = ggplot2::aes(x = lambda, y = value, group = name)) +
      ggplot2::geom_line(colour = "#008080")+
      ggplot2::ggtitle("Regularized Parameters")
    
    fitPlot <- ggplot2::ggplot(data = fits,
                               mapping = ggplot2::aes(x = lambda, y = cvfit), title = "cv-fit") +
      ggplot2::geom_line()
    
  }else{
    fits <- fits[fits$alpha == alpha,, drop = FALSE]
    parameters <- cbind(
      lambda = parameters$lambda,
      alpha = parameters$alpha,
      coef(x, alpha = alpha)
    )
    parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@parameterLabels)
    parametersLong$regularized <- parametersLong$name %in% x@regularized
    
    parameterPlot <- ggplot2::ggplot(data = parametersLong,
                    mapping = ggplot2::aes(x = lambda, 
                                           y = value, 
                                           group = name, 
                                           color = regularized)) +
      ggplot2::geom_line()+ 
      ggplot2::scale_color_manual(values=c("FALSE"="gray","TRUE"="#008080")) +
      ggplot2::ggtitle("Parameter Estimates")
    
    fitPlot <- ggplot2::ggplot(data = fits,
                               mapping = ggplot2::aes(x = lambda, y = cvfit), title = "cv-fit") +
      ggplot2::geom_line()
  }
  
  legendIs <- cowplot::get_legend(parameterPlot)
  
  plotIs <- cowplot::plot_grid(
    parameterPlot + ggplot2::theme(legend.position="none"),
    fitPlot,
    align = 'vh',
    hjust = -1,
    nrow = 1
  )
  cowplot::plot_grid(plotIs, legendIs, rel_widths = c(3, .4))
})