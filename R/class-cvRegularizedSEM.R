#' Class for cross-validated regularized SEM
#' @slot parameters data.frame with parameter estimates for the best combination of the
#' tuning parameters
#' @slot transformations transformed parameters
#' @slot cvfits data.frame with all combinations of the
#' tuning parameters and the sum of the cross-validation fits
#' @slot parameterLabels character vector with names of all parameters
#' @slot regularized character vector with names of regularized parameters
#' @slot cvfitsDetails data.frame with cross-validation fits for each subset
#' @slot subsets matrix indicating which person is in which subset
#' @slot subsetParameters optional: data.frame with parameter estimates for all
#' combinations of the tuning parameters in all subsets
#' @slot misc list with additional return elements
setClass(Class = "cvRegularizedSEM",
         representation = representation(
           parameters="data.frame",
           transformations="data.frame",
           cvfits = "data.frame",
           parameterLabels = "character",
           regularized = "character",
           cvfitsDetails="data.frame", 
           subsets = "matrix",
           subsetParameters = "data.frame",
           misc = "list"
         )
)

#' Show method for objects of class \code{cvRegularizedSEM}.
#' 
#' @param object object of class cvRegularizedSEM
#' @return No return value, just prints estimates
#' @export
setMethod("show", "cvRegularizedSEM", function (object) {
  bestFit <- unlist(object@parameters)
  tuningParameters <- bestFit[!names(bestFit) %in% object@parameterLabels]
  cat(paste0("Best parameter estimates observed for: ",
             paste0(names(tuningParameters), " = ", tuningParameters, collapse = "; "), 
             ":\n"))
  print(bestFit[object@parameterLabels])
})

#' summary method for objects of class \code{cvRegularizedSEM}.
#' 
#' @param object object of class cvRegularizedSEM
#' @return No return value, just prints estimates
#' @export
setMethod("summary", "cvRegularizedSEM", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Exact Cross Validation Results ####\n\n"))
  cat("regularized parameters: ")
  cat(paste0(object@regularized, collapse = ", "))
  cat("\n\n")
  cat("Final parameter estimates based on best out-of-sample fit:\n")
  print(coef(object))
  cat("\n\n")
  cat(paste0("- Use coef(", modelName, 
             ") to get the final parameter estimates of the model.\n\n"))
  cat(paste0("- Use plot(", modelName, 
             ") to plot the cross-validation fit of the model.\n\n"))
  cat(paste0("- Cross-validation fits can be accessed with ", modelName, 
             "@cvfits.\n\n"))
  cat("################################################\n")
})

#' coef
#' 
#' Returns the parameter estimates of an cvRegularizedSEM
#'  
#' @param object object of class cvRegularizedSEM
#' @param ... not used
#' @return the parameter estimates of an cvRegularizedSEM
#' @importMethodsFrom stats coef
#' @export
setMethod("coef", "cvRegularizedSEM", function (object, ...) {
  
  tuningParameters <- object@parameters[, !colnames(object@parameters) %in% object@parameterLabels,drop=FALSE] 
  estimates <- as.matrix(object@parameters[,object@parameterLabels,drop=FALSE])
  
  coefs <- new("lessSEMCoef")
  coefs@tuningParameters <- tuningParameters
  coefs@estimates <- estimates
  
  return(coefs)
})

#' plots the cross-validation fits
#' 
#' @param x object of class cvRegularizedSEM
#' @param y not used
#' @param ... not used
#' @importMethodsFrom graphics plot
#' @return either an object of ggplot2 or of plotly
#' @export
setMethod("plot", 
          signature = c(x = "cvRegularizedSEM",
                        y = "missing"), 
          definition = function (x, y, ...) 
          {
            
            fits <- x@cvfits
            tuningParameters <- fits[,colnames(fits)!= "cvfit",drop=FALSE]
            tuningParameters <- tuningParameters[,apply(tuningParameters,2,function(x) length(unique(x)) > 1),drop=FALSE]
            
            nTuning <- ncol(tuningParameters)
            
            if(nTuning > 2) 
              stop("Plotting currently only supported for up to 2 tuning parameters")
            if(nTuning == 2 & !("plotly" %in% rownames(utils::installed.packages())))
              stop("Plotting more than one tuning parameter requires the package plotly")
            
            if(nTuning == 1){
              # .data[[v1]], .data[["y"]])
              return(
                ggplot2::ggplot(data = fits,
                                mapping = ggplot2::aes(
                                  x = .data[[colnames(tuningParameters)]], 
                                  y = .data[["cvfit"]])) +
                  ggplot2::geom_line(colour = "#008080") +
                  ggplot2::ggtitle("Regularized Parameters")
              )
              
            }else{
              fits$tp1 <- unlist(fits[,colnames(tuningParameters)[1]])
              fits$tp2 <- unlist(fits[,colnames(tuningParameters)[2]])
              return(
                plotly::layout(
                  plotly::plot_ly(fits, 
                                  x = ~tp1, y = ~tp2, z = ~cvfit, 
                                  type = 'scatter3d', mode = 'lines',
                                  opacity = 1,
                                  split = ~tp2,
                                  line = list(width = 6, 
                                              reverscale = FALSE)), 
                  scene = list(xaxis = list(title = colnames(tuningParameters)[1]), 
                               yaxis = list(title = colnames(tuningParameters)[2]))
                )
              )
              
            }
          })