#' Class for cross-validated regularized SEM
#' @slot parameters data.frame with parameter estimates for the best combination of the
#' tuning parameters
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
#' @docType methods
#' @name show-cvRegularizedSEM
#' @rdname show-cvRegularizedSEM
#' @aliases show-cvRegularizedSEM show,cvRegularizedSEM-method
#' 
#' @param object object of class cvRegularizedSEM
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
#' @docType methods
#' @name summary-cvRegularizedSEM
#' @rdname summary-cvRegularizedSEM
#' @aliases summary-cvRegularizedSEM summary,cvRegularizedSEM-method
#' 
#' @param object object of class cvRegularizedSEM
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
#' @returns the parameter estimates of an cvRegularizedSEM
#' @export
setMethod("coef", "cvRegularizedSEM", function (object) {
  return(unlist(object@parameters[,object@parameterLabels]))
})

#' plots the regularized and unregularized parameters as well as the cross-validation fits for all levels of lambda
#'
#' @docType methods
#' @name plot-cvRegularizedSEM
#' @rdname plot-cvRegularizedSEM
#' @aliases plot-cvRegularizedSEM plot,cvRegularizedSEM-method
#' 
#' 
#' @param x object of class cvRegularizedSEM
setMethod("plot", "cvRegularizedSEM", function (x) {
  
  fits <- x@cvfits
  tuningParameters <- fits[,colnames(fits)!= "cvfit",drop=FALSE]
  tuningParameters <- tuningParameters[,apply(tuningParameters,2,function(x) length(unique(x)) > 1),drop=FALSE]
  
  nTuning <- ncol(tuningParameters)
  
  if(nTuning > 2) 
    stop("Plotting currently only supported for up to 2 tuning parameters")
  if(nTuning == 2 & !("plotly" %in% rownames(utils::installed.packages())))
    stop("Plotting more than one tuning parameter requires the package plotly")
  
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
})