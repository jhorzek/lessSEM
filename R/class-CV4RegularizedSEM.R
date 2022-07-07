setClass(Class = "CV4RegularizedSEM",
         representation = representation(
           parameters="data.frame",
           cvfits = "data.frame",
           parameterLabels = "character",
           regularized = "character",
           cvfitsDetails="data.frame", 
           subsets = "matrix",
           subsetParameters = "data.frame"
         )
)

#' show
#' 
#' @export
setMethod("show", "CV4RegularizedSEM", function (object) {
  bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
  tuningParameters <- object@parameters[,!colnames(object@parameters)%in%object@parameterLabels]
  paste0(colnames(tuningParameters), " = ", tuningParameters[bestFit,])
  cat(paste0("Best parameter estimates observed for: ",
             paste0(colnames(tuningParameters), " = ", tuningParameters[bestFit,], collapse = "; "), 
             ":\n"))
  
  parameters <- unlist(object@parameters[bestFit,object@parameterLabels])
  print(parameters)
})

#' summary
#' 
#' @export
setMethod("summary", "CV4RegularizedSEM", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Exact Cross Validation Results ####\n\n"))
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
#' Returns the parameter estimates of an CV4RegularizedSEM
#'  
#' @param object object of class regularizedSEM
#' @param lambda numeric: value of lambda for which parameters should be returned. Must be present in object
#' @param alpha numeric: value of alpha for which parameters should be returned. Must be present in object
#' @param rule one of "min" or "1sd". With "min", the parameters of the model with the lowest CV-fit value will be returned.
#' With "1sd", the 1 standard deviation rule is used and the parameters with the larges lambda and a fit within 1 sd of the minimum will be returned.
#' 
#' @export
setMethod("coef", "CV4RegularizedSEM", function (object, rule = "min") {
  
    if(!rule %in% c("min", "1sd")) stop("rule must be one of: 'min', '1sd'.")
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
    
  return(pars)
})

#' plot
#' 
#' plots the regularized and unregularized parameters as well as the cross-validation fits for all levels of lambda
#' 
#' @param x object of class regularizedSEM
#' @param what either "parameters" or "fit"
#' @param regularizedOnly boolean: should only regularized parameters be plotted?``
#' @export
setMethod("plot", "CV4RegularizedSEM", function (x, what = "parameters", regularizedOnly = TRUE) {
  
  parameters <- x@parameters
  fits <- x@cvfits
  tuningParameters <- x@parameters[,!colnames(x@parameters)%in%x@parameterLabels,drop=FALSE]
  tuningParameters <- tuningParameters[,apply(tuningParameters,2,function(x) length(unique(x)) > 1),drop=FALSE]
  
  nTuning <- ncol(tuningParameters)
  
  if(nTuning > 2) 
    stop("Plotting currently only supported for up to 2 tuning parameters")
  if(nTuning == 2 & !("plotly" %in% rownames(installed.packages())))
    stop("Plotting more than one tuning parameter requires the package plotly")
  
  if(what == "parameters"){
    
    if(regularizedOnly){
      
      parameters <- cbind(
        tuningParameters,
        parameters[,x@regularized, drop = FALSE]
      )
      parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@regularized)
      
    }else{
      
      parameters <- cbind(
        tuningParameters,
        parameters
      )
      parametersLong <- tidyr::pivot_longer(data = parameters, cols = x@parameterLabels)
      
    }
    
    if(nTuning == 1){
      
      ggplot2::ggplot(data = parametersLong,
                      mapping = ggplot2::aes_string(x = colnames(tuningParameters), 
                                                    y = "value", 
                                                    group = "name")) +
        ggplot2::geom_line(colour = "#008080")+
        ggplot2::ggtitle("Regularized Parameters")
      
    }else{
      parametersLong$name <- paste0(parametersLong$name, 
                                    "_", 
                                    unlist(parametersLong[,colnames(tuningParameters)[2]]))
      parametersLong$tp1 <- unlist(parametersLong[,colnames(tuningParameters)[1]])
      parametersLong$tp2 <- unlist(parametersLong[,colnames(tuningParameters)[2]])
      plt <- plotly::plot_ly(parametersLong, 
                             x = ~tp1, y = ~tp2, z = ~value, 
                             type = 'scatter3d',
                             mode = 'lines',
                             opacity = 1,
                             color = ~name,
                             split = ~tp2,
                             line = list(width = 6, 
                                         reverscale = FALSE)
      )
      print(plt)
      
    }
    
  }else{
    if(nTuning == 1){
      
      ggplot2::ggplot(data = fits,
                      mapping = ggplot2::aes_string(x = colnames(tuningParameters), 
                                                    y = "cvfit")) +
        ggplot2::geom_line(colour = "#008080")+
        ggplot2::ggtitle("Regularized Parameters")
      
    }else{
      fits$tp1 <- unlist(fits[,colnames(tuningParameters)[1]])
      fits$tp2 <- unlist(fits[,colnames(tuningParameters)[2]])
      plotly::plot_ly(fits, 
                      x = ~tp1, y = ~tp2, z = ~cvfit, 
                      type = 'scatter3d', mode = 'lines',
                      opacity = 1,
                      split = ~tp2,
                      line = list(width = 6, 
                                  reverscale = FALSE))
      
    }
  }
})