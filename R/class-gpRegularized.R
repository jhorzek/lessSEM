#' Class for regularized model using general purpose optimization interface
#' @slot penalty penalty used (e.g., "lasso")
#' @slot parameters data.frame with all parameter estimates
#' @slot fits data.frame with all fit results
#' @slot parameterLabels character vector with names of all parameters
#' @slot weights vector with weights given to each of the parameters in the penalty
#' @slot regularized character vector with names of regularized parameters
#' @slot internalOptimization list of elements used internally
#' @slot inputArguments list with elements passed by the user to the general
#' purpose optimizer
#' @export
setClass(Class = "gpRegularized",
         representation = representation(
           penalty = "character",
           parameters="data.frame",
           fits="data.frame", 
           parameterLabels = "character",
           weights = "numeric",
           regularized = "character",
           internalOptimization="list", 
           inputArguments="list"
         )
)

#' show
#' 
#' @param object object of class gpRegularized
#' @return No return value, just prints estimates
#' @export
setMethod("show", "gpRegularized", function (object) {
  #modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class gpRegularized with ",object@penalty, " penalty ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat(paste0("- Use coef(object) to get the parameter estimates of the model. With coef(object, lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Use plot(object) to plot the parameter estimates of the model.\n\n"))
  cat(paste0("- Information criteria can be computed with AIC(object) or BIC(object).\n\n"))
  cat("################################################\n")
})

#' summary
#' 
#' @param object object of class gpRegularized
#' @param ... not used
#' @return No return value, just prints estimates
#' @export
setMethod("summary", "gpRegularized", function (object, ...) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Model of class gpRegularized with ",object@inputArguments$penalty, " penalty ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat(paste0("- Use coef(", modelName, 
             ") to get the parameter estimates of the model. With coef(", 
             modelName, "lambda = x, delta = y) parameters estimates at the values x and y for lambda and delta can be extracted.\n\n"))
  cat(paste0("- Use plot(", modelName, 
             ") to plot the parameter estimates of the model.\n\n"))
  cat(paste0("- Information criteria can be computed with AIC(", modelName, 
             ") or BIC(", modelName, 
             ").\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of a gpRegularized
#' 
#' @param object object of class gpRegularized
#' @param ... criterion can be one of: "AIC", "BIC". If set to NULL, all parameters will be returned
#' @return parameter estimates
#' @export
setMethod("coef", "gpRegularized", function (object, ...) {
  dotdotdot <- list(...)
  if("criterion" %in% names(dotdotdot)){
    criterion <- dotdotdot$criterion
  }else{
    criterion <- NULL
  }
  
  tuningParameters <- object@parameters[, !colnames(object@parameters) %in% object@parameterLabels,drop=FALSE] 
  estimates <- as.matrix(object@parameters[,object@parameterLabels,drop=FALSE])
  
  if(!is.null(criterion) && criterion %in% c("AIC", "BIC")){
    if(length(unique(object@fits$nonZeroParameters)) == 1) 
      stop("Selection by criterion currently only supported for sparsity inducing penalties. Either none of your parameters was zeroed or the penalty used does not induce sparsity.")
    if(criterion == "AIC"){
      AICs <- AIC(object)
      bestAIC <- which(AICs$AIC == min(AICs$AIC))[1]
      
      coefs <- new("lessSEMCoef")
      coefs@tuningParameters <- tuningParameters[bestAIC,,drop = FALSE]
      coefs@estimates <- estimates[bestAIC,,drop = FALSE]
      coefs@transformations <- matrix(nrow = 0, ncol = 0)
      
      return(coefs) 
    }
    
    if(criterion == "BIC"){
      BICs <- BIC(object)
      bestBIC <- which(BICs$BIC == min(BICs$BIC))[1]
      
      coefs <- new("lessSEMCoef")
      coefs@tuningParameters <- tuningParameters[bestBIC,,drop = FALSE]
      coefs@estimates <- estimates[bestBIC,,drop = FALSE]
      coefs@transformations <- matrix(nrow = 0, ncol = 0)
      
      return(coefs)
    }
    
  }
  
  coefs <- new("lessSEMCoef")
  coefs@tuningParameters <- tuningParameters
  coefs@estimates <- estimates
  coefs@transformations <- matrix(nrow = 0, ncol = 0)
  
  return(coefs)
})

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class gpRegularized
#' @param ... not used
#' @param k multiplier for number of parameters
#' @return data frame with fit values, appended with AIC
#' @export
setMethod("AIC", "gpRegularized", function (object, ..., k = 2) {
  if(!object@penalty %in% c("lasso", "adaptiveLasso", "cappedL1", "mcp", "scad"))
    stop("AIC not supported for this penalty.")
  
  fits <- object@fits
  fits$AIC <- fits$m2LL + k*fits$nonZeroParameters
  
  return(fits)
})

#' BIC
#' 
#' returns the BIC
#' 
#' @param object object of class gpRegularized
#' @param ... not used
#' @return data frame with fit values, appended with BIC
#' @export
setMethod("BIC", "gpRegularized", function (object, ...) {
  N <- nrow(lavaan::lavInspect(object@inputArguments$lavaanModel, "data"))
  fits <- object@fits
  
  if(!object@penalty %in% c("lasso", "adaptiveLasso", "cappedL1", "mcp", "scad"))
    stop("BIC not supported for this penalty.")
  
  fits <- object@fits
  fits$BIC <- fits$m2LL + log(N)*fits$nonZeroParameters
  
  return(fits)
  
})

#' plots the regularized and unregularized parameters for all levels of lambda
#' 
#' @param x object of class gpRegularized
#' @param y not used
#' @param ... use regularizedOnly=FALSE to plot all parameters
#' @return either an object of ggplot2 or of plotly
#' @export
setMethod("plot",
          c(x = "gpRegularized", y = "missing"), 
          function (x, y, ...)
          {
            if("regularizedOnly" %in% names(list(...))){
              regularizedOnly <- list(...)$regularizedOnly
            }else{
              regularizedOnly <- TRUE
            }
            parameters <- x@parameters
            tuningParameters <- x@parameters[,!colnames(x@parameters)%in%x@parameterLabels,drop=FALSE]
            tuningParameters <- tuningParameters[,apply(tuningParameters,2,function(x) length(unique(x)) > 1),drop=FALSE]
            
            nTuning <- ncol(tuningParameters)
            
            if(nTuning > 2) 
              stop("Plotting currently only supported for up to 2 tuning parameters")
            if(nTuning == 2 & !("plotly" %in% rownames(utils::installed.packages())))
              stop("Plotting more than one tuning parameter requires the package plotly")
            
            
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
              
              return(
                ggplot2::ggplot(data = parametersLong,
                                mapping = ggplot2::aes(
                                  x = .data[[colnames(tuningParameters)]], 
                                  y = .data[["value"]],
                                  group = .data[["name"]])) +
                  ggplot2::geom_line(colour = "#008080")+
                  ggplot2::ggtitle("Regularized Parameters")
              )
              
            }else{
              parametersLong$name <- paste0(parametersLong$name, 
                                            "_", 
                                            unlist(parametersLong[,colnames(tuningParameters)[2]]))
              parametersLong$tp1 <- unlist(parametersLong[,colnames(tuningParameters)[1]])
              parametersLong$tp2 <- unlist(parametersLong[,colnames(tuningParameters)[2]])
              plt <- plotly::layout(
                plotly::plot_ly(parametersLong, 
                                x = ~tp1, y = ~tp2, z = ~value, 
                                type = 'scatter3d',
                                mode = 'lines',
                                opacity = 1,
                                color = ~name,
                                split = ~tp2,
                                line = list(width = 6, 
                                            reverscale = FALSE)
                ), 
                scene = list(xaxis = list(title = colnames(tuningParameters)[1]), 
                             yaxis = list(title = colnames(tuningParameters)[2]))
              )
              return(plt)
              
            }
            
          })
