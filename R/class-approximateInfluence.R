setClass("approximateInfluence",
         representation = representation(
           subsets="matrix",
           tuningParameters = "data.frame",
           subsetParameters="data.frame",
           subsetFits = "data.frame",
           parameters = "data.frame",
           parameterLabels = "character",
           regularized  = "character"
         ))

#' plot
#' 
#' 
#' @param x object of class approximateInfluence
#' @param alpha the alpha tuning parameter value for which the data should be plotted
#' @param regularizedOnly boolean: Should only the regularized parameters be plotted?
#' @param interactive boolean: If set to TRUE, plotly will be used to make the plot interactive. This can provide more information into which person may influence the results substantially 
#' @export
setMethod("plot", "approximateInfluence", function (x, alpha = NULL, regularizedOnly = TRUE, interactive = FALSE) {
  parameters <- x@parameters
  if(is.null(alpha)) {
    alpha <- unique(parameters$alpha)[1]
    
    if(length(unique(parameters$alpha)) != 1) warning(paste0("Models for different values of alpha were fitted, but not alpha was provided to the plotting function. Plotting for alpha = ", 
                                                             alpha, ". You can specify specific alpha values with, for instance, plot(approximateInfluence, alpha = .02)"))
  }
  if(!alpha %in% parameters$alpha){
    stop(paste0("alpha = ", alpha, " is was not tested in the model. Only found the following alpha values: "),
         paste0(unique(parameters$alpha), collapse= ", "), "."
    )
  }
  
  if(regularizedOnly){
    selectParameters <- x@regularized
  }else{
    selectParameters <- x@parameterLabels
  }
  
  parameters <- parameters[parameters$alpha == alpha,]
  parameters <- parameters[,c("lambda", selectParameters), drop = FALSE]
  
  subsetParameters <- x@subsetParameters[x@subsetParameters$alpha == alpha,]
  subsetParameters <- subsetParameters[,c("removedSubset", "lambda", selectParameters), drop = FALSE]
  
  parametersLong <- tidyr::pivot_longer(data = parameters, cols = selectParameters)
  subsetParametersLong <- tidyr::pivot_longer(data = subsetParameters, cols = selectParameters)
  subsetParametersLong$removedSubset <- as.factor(subsetParametersLong$removedSubset)
  
  plt <- ggplot2::ggplot(parametersLong, 
                         mapping = ggplot2::aes(x = lambda, y = value, color = name)) +
    ggplot2::geom_line() + 
    ggplot2::ggtitle("Regularized Parameters")
  
  plt <- plt + 
    ggplot2::geom_point(data = subsetParametersLong, ggplot2::aes(x = lambda, y = value, color = name, group = removedSubset))
  
  if(interactive){
    if(!"plotly" %in% rownames(installed.packages())){
      warning("Interactive plots require plotly to be installed (install.packages('plotly')). Returning non-interactive plot.")
      print(plt)
    }else{
      plotly::ggplotly(plt)
    }
  }else{
    print(plt)
  }
})

setClass("approximateInfluenceFitIndices",
         representation = representation(
           informationCriterion = "character",
           values = "data.frame",
           subsetParameters="data.frame"
         ))

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class approximateInfluence
#' @export
setMethod("AIC", "approximateInfluence", function (object) {
  AICs <- object@subsetFits
  AICs$nonZeroParameters <- apply(object@subsetParameters[,object@parameterLabels] != 0,1,sum)
  AICs$value <- AICs$fit + 2*AICs$nonZeroParameters
  return(new("approximateInfluenceFitIndices",
             informationCriterion = "AIC",
             values = AICs,
             subsetParameters = object@subsetParameters
  ))
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class approximateInfluence
#' @export
setMethod("BIC", "approximateInfluence", function (object) {
  BICs <- object@subsetFits
  BICs$nonZeroParameters <- apply(object@subsetParameters[,object@parameterLabels] != 0,1,sum)
  BICs$value <- BICs$fit + log(apply(!object@subsets[object@subsetParameters$removedSubset,],1,sum))*BICs$nonZeroParameters
  return(new("approximateInfluenceFitIndices",
             informationCriterion = "BIC",
             values = BICs,
             subsetParameters = object@subsetParameters
  ))
})

#' coef
#' 
#' Returns the parameter estimates for selected subsets. By default, the best parameter estimates for each subset are returned.
#' 
#' @param object object of class approximateInfluenceFitIndices
#' @param lambda for which lambda should the parameters be returned? If set to NULL, lambda is selected based on the information criterion
#' @param removedSubset for which removed subset should the parameters be returned? If set to NULL, parameters are returned for all removed subsets.
#' @export
setMethod("coef", "approximateInfluenceFitIndices", function (object, lambda = NULL, removedSubset = NULL) {
  alpha <- unique(object@values$alpha)[1]
  if(length(unique(object@values$alpha)) != 1 || alpha != 1) stop("Currently only supported for lasso type penalties.")
  
  values <- object@values
  parameters <- object@subsetParameters
  if(!is.null(removedSubset)){
    values <- values[values$removedSubset %in% removedSubset,]
    parameters <- parameters[parameters$removedSubset %in% removedSubset,]
  }
  
  if(!is.null(lambda)){
    if(!lambda %in% unique(values$lambda)) stop("Could not find requested lambda in the information criteria. Please check object@values$lambda.")
    values <- values[values$lambda %in% lambda,]
    parameters <- parameters[parameters$lambda %in% lambda,]
    return(list(
      "fitIndices" = values,
      "parameters" = parameters
    ))
  }
  
  # get minimum of information criterion for each person
  whichIC <- c(by(values$value, list(values$removedSubset), FUN = function(z) which(z == min(z))[1]))
  
  finalParameters <- matrix(NA, nrow = length(whichIC), ncol = ncol(parameters),
                            dimnames = list(NULL, colnames(parameters)))
  for(i in values$removedSubset){
    parameters_i <- parameters[parameters$removedSubset == i,]
    finalParameters[i,] <- unlist(parameters_i[whichIC[i],colnames(finalParameters)])
  }
  
  return(finalParameters)
})

#' mean
#' 
#' Returns the mean information criterion value over all subsets
#' 
#' @param x object of class approximateInfluenceFitIndices
#' @export
setMethod("mean", "approximateInfluenceFitIndices", function (x) {
  alpha <- unique(x@values$alpha)[1]
  if(length(unique(x@values$alpha)) != 1 || alpha != 1) stop("Currently only supported for lasso type penalties.")
  
  values <- c(by(x@values$value, list(x@values$lambda), mean))
  ICMean <- data.frame(
    lambda = unique(x@values$lambda),
    value = values
  )
  attr(ICMean, "informationCriterion") <- x@informationCriterion
  return(ICMean)
})

#' median
#' 
#' Returns the median information criterion value over all subsets
#' 
#' @param x object of class approximateInfluenceFitIndices
#' @export
setMethod("median", "approximateInfluenceFitIndices", function (x, na.rm = NULL) {
  alpha <- unique(x@values$alpha)[1]
  if(length(unique(x@values$alpha)) != 1 || alpha != 1) stop("Currently only supported for lasso type penalties.")
  
  values <- c(by(x@values$value, list(x@values$lambda), median))
  ICMean <- data.frame(
    lambda = unique(x@values$lambda),
    value = values
  )
  attr(ICMean, "informationCriterion") <- x@informationCriterion
  return(ICMean)
})

#' plot
#' 
#' @param x object of class approximateInfluenceFitIndices
#' @param alpha the alpha tuning parameter value for which the fit index should be plotted
#' @param interactive boolean: If set to TRUE, plotly will be used to make the plot interactive. This can provide more information into which person may influence the results substantially 
#' @export
setMethod("plot", "approximateInfluenceFitIndices", function (x, alpha = NULL, interactive = FALSE) {
  
  if(is.null(alpha)) {
    alpha <- unique(x@values$alpha)[1]
    
    if(length(unique(x@values$alpha)) != 1) warning(paste0("Models for different values of alpha were fitted, but not alpha was provided to the plotting function. Plotting for alpha = ", 
                                                           alpha, ". You can specify specific alpha values with, for instance, plot(approximateInfluence, alpha = .02)"))
  }
  if(!alpha %in% x@values$alpha){
    stop(paste0("alpha = ", alpha, " is was not tested in the model. Only found the following alpha values: "),
         paste0(unique(x@values$alpha), collapse= ", "), "."
    )
  }
  
  data2Plot <- data.frame(
    value = x@values$value[x@values$alpha == alpha],
    lambda = x@values$lambda[x@values$alpha == alpha],
    removedSubset = as.factor(x@values$removedSubset[x@values$alpha == alpha])
  )
  pltTitle <- paste0("Influence on ", x@informationCriterion)
  
  plt <- ggplot2::ggplot(data2Plot, 
                         mapping = ggplot2::aes(x = lambda, y = value, color = removedSubset)) +
    ggplot2::geom_line() + 
    ggplot2::ggtitle(pltTitle) + 
    ggplot2::theme(legend.position = "none")
  
  if(interactive){
    if(!"plotly" %in% rownames(installed.packages())){
      warning("Interactive plots require plotly to be installed (install.packages('plotly')). Returning non-interactive plot.")
      print(plt)
    }else{
      plotly::ggplotly(plt)
    }
  }else{
    print(plt)
  }
  
})