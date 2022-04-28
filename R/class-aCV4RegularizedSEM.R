setClass(Class = "aCV4RegularizedSEM",
         representation = representation(
           parameters="data.frame",
           cvfits = "data.frame",
           parameterLabels = "character",
           regularized = "character",
           cvfitsDetails="data.frame", 
           subsets = "matrix"
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
#' # References
#' 
#' Ternès, N., Rotolo, F., & Michiels, S. (2016). Empirical extensions of the lasso penalty to reduce the false discovery rate in high-dimensional Cox regression models: Extended lasso regression reducing false positives in Cox models. Statistics in Medicine, 35(15), 2561–2573. https://doi.org/10.1002/sim.6927
#' 
#' @param object object of class regularizedSEM
#' @param lambda numeric: value of lambda for which parameters should be returned. Must be present in object
#' @param alpha numeric: value of alpha for which parameters should be returned. Must be present in object
#' @param rule one of "min", "1sd", and "penalized". With "min", the parameters of the model with the lowest CV-fit value will be returned.
#' With "1sd", the 1 standard deviation rule is used and the parameters with the larges lambda and a fit within 1 sd of the minimum will be returned.
#' With "penalized", parameters with the additional reward for sparsity proposed by (Ternès et al., 2016) will be returned
#' 
#' @export
setMethod("coef", "aCV4RegularizedSEM", function (object, lambda = NULL, alpha = NULL, rule = "min") {
  
  if(is.null(lambda) && is.null(alpha)) {
    # using rule
    if(!rule %in% c("min", "1sd", "penalized")) stop("rule must be one of: 'min', '1sd', 'penalized'.")
    if(rule == "min"){
      # return parameters of model with best cv fit
      bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
      parameters <- unlist(object@parameters[bestFit,object@parameterLabels])
      return(parameters)
    }
    if(rule == "1sd"){
      # return parameters of model with highest lambda within 1 sd of minimum
      oneSd <- sd(object@cvfits$cvfit)
      withinOneSd <- which(abs(object@cvfits$cvfit - min(object@cvfits$cvfit)) < oneSd)
      select <- max(withinOneSd)
      parameters <- unlist(object@parameters[select,object@parameterLabels])
      return(parameters)
    }
    if(rule == "penalized"){
      # Ternès, N., Rotolo, F., & Michiels, S. (2016). Empirical extensions of the lasso penalty to reduce the false discovery rate in high-dimensional Cox regression models: Extended lasso regression reducing false positives in Cox models. Statistics in Medicine, 35(15), 2561–2573. https://doi.org/10.1002/sim.6927
      if(any(object@cvfits$alpha != 1)) stop("rule = 'penalized' only supported for lasso and adaptive lasso")
      
      cvFits <- object@cvfits$cvfit
      parameters <- object@parameters[,object@regularized]
      parameterWasZeroed <- apply(parameters==0,2,sum) > 0 # check if all parameters have been zeroed for at least one lambda
      nNotZeroed <- apply(parameters != 0, 1, sum)
      
      if(!all(parameterWasZeroed)) {
        warning("rule = 'penalized' requires that ALL regularized parameters have been set to zero for at least one lambda. This is not satisfied in the object passed here. Using the maximal lambda instead.")
        lambdaMax <- which(object@cvfits$lambda == max(object@cvfits$lambda))
      }else{
        # the lambda for which all parameters were first zeroed provides the upper bound
        lambdaMax <- min(which(nNotZeroed == 0))
      }
      lambdaMin <- which.min(cvFits) # the lambda which has the best cv fit provides the lower bound
      
      bestCVFit <- min(cvFits)
      sparseCVFit <- cvFits[lambdaMax]
      
      bestNNonzeroParameter <- nNotZeroed[lambdaMin] 
      
      penalty <- ((sparseCVFit - bestCVFit)/bestNNonzeroParameter)*nNotZeroed[lambdaMin:lambdaMax] # note: we use the -2 log-Likelihood; our penalty term is therefore different from
      # Ternès et al. (2016)
      
      cvFitsPenalized <- cvFits
      cvFitsPenalized[] <- Inf
      cvFitsPenalized[lambdaMin:lambdaMax] <- cvFits[lambdaMin:lambdaMax] + penalty
      bestFit <- which.min(cvFitsPenalized)
      
      parameters <- unlist(object@parameters[bestFit,object@parameterLabels])
      return(parameters)
    }
   
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