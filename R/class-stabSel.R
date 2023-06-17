#' Class for stability selection
#' @slot regularized names of regularized parameters
#' @slot tuningParameters data.frame with tuning parameter values
#' @slot stabilityPaths matrix with percentage of parameters being non-zero
#' averaged over all subsets for each setting of the tuning parameters
#' @slot percentSelected percentage with which a parameter was selected over all
#' tuning parameter settings
#' @slot selectedParameters final selected parameters
#' @slot settings internal
#' @export
setClass("stabSel",
         representation = representation(
           regularized = "character",
           tuningParameters = "data.frame",
           stabilityPaths = "matrix",
           percentSelected = "numeric",
           selectedParameters = "logical",
           settings = "list"
         ))

#' show
#' 
#' @param object object of class stabSel
#' @return No return value, just prints estimates
setMethod("show", "stabSel", function (object) {
  cat("Results for stability selection\n")
  cat(paste0(rep("-", nchar("Threshold for adding parameter:   100 %")), collapse = ""), "\n")
  cat("Number of subsamples:             ", object@settings$numberOfSubsamples, "\n")
  cat("Sample size in subsamples:        ", object@settings$subsampleSize, "\n") 
  cat("Threshold for adding parameter:   ", object@settings$threshold, "%", "\n") 
  cat(paste0(rep("-", nchar("Threshold for adding parameter:   100 %")), collapse = ""), "\n")
  cat("\n\nUnregularized parameters:\n")
  cat(paste0(names(object@selectedParameters)[!names(object@selectedParameters) %in% object@regularized],
      collapse = ", "))
  
  cat("\n\nRegularized parameters above threshold:\n")
  cat(paste0(names(object@selectedParameters)[(names(object@selectedParameters) %in% object@regularized) &
                                                object@selectedParameters],
             collapse = ", "))
  
  cat("\n\nRemoved parameters:\n")
  cat(paste0(names(object@selectedParameters)[(names(object@selectedParameters) %in% object@regularized) &
                                                !object@selectedParameters],
             collapse = ", "))
  cat("\n")
})

#' plots the regularized and unregularized parameters for all levels of the tuning parameters
#' 
#' @param x object of class stabSel
#' @param y not used
#' @param ... use regularizedOnly=FALSE to plot all parameters
#' @return either an object of ggplot2 or of plotly
#' @export
setMethod("plot", 
          c(x = "stabSel", y = "missing"), 
          function (x, y, ...) {
            if("regularizedOnly" %in% names(list(...))){
              regularizedOnly <- list(...)$regularizedOnly
            }else{
              regularizedOnly <- TRUE
            }
            
            parameters <- x@stabilityPaths
            tuningParameters <- x@tuningParameters
            tuningParameters <- tuningParameters[,apply(tuningParameters,2,function(x) length(unique(x)) > 1),drop=FALSE]
            
            nTuning <- ncol(tuningParameters)
            
            if(nTuning > 2) 
              stop("Plotting currently only supported for up to 2 tuning parameters")
            if(nTuning == 2 & !requireNamespace("plotly", quietly = TRUE))
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
              parametersLong <- tidyr::pivot_longer(data = parameters, cols = colnames(parameters))
              
            }
            
            if(nTuning == 1){
              
              return(
                ggplot2::ggplot(data = parametersLong,
                                mapping = ggplot2::aes(
                                  x = .data[[colnames(tuningParameters)]], 
                                  y = .data[["value"]],
                                  group = .data[["name"]])) +
                  ggplot2::ylab("% of models with this parameter being non-zero") +
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
