# Note: we define a custom logLik - Function because the generic one is 
# using df = number of parameters which might be confusing.
setClass("approximateInfluence",
         representation = representation(
           subsets="matrix",
           tuningParameters = "data.frame",
           subsetParameters="data.frame",
           subsetFits = "data.frame",
           parameterLabels = "character",
           regularized  = "character"
         ))

#' AIC
#' 
#' returns the AIC
#' 
#' @param object object of class approximateInfluence
#' @export
setMethod("AIC", "approximateInfluence", function (object) {
  AICs <- object@subsetFits
  AICs$nonZeroParameters <- apply(object@subsetParameters[,object@parameterLabels] == 0,1,sum)
  AICs$AIC <- AICs$fit + 2*AICs$nonZeroParameters
  return(AICs)
})

#' AIC
#' 
#' returns the BIC
#' 
#' @param object object of class approximateInfluence
#' @export
setMethod("BIC", "approximateInfluence", function (object) {
  BIC <- object@subsetFits
  BIC$nonZeroParameters <- apply(object@subsetParameters[,object@parameterLabels] == 0,1,sum)
  BIC$BIC <- BIC$fit + log(nrow(object@subsets))*BIC$nonZeroParameters
  return(BIC)
})

#' plot
#' 
#' 
#' @param object object of class approximateInfluence
#' @export
setMethod("plot", "approximateInfluence", function (x, alpha = NULL) {
  if(is.null(alpha)) alpha <- 1
  parameters <- x@subsetParameters[x@subsetParameters$alpha == alpha,]
  parameters <- parameters[,c("removedSubset", "lambda", x@regularized), drop = FALSE]
  
  parametersDescr <- c()
  for(l in unique(parameters$lambda)){
    parametersDescr <- rbind(parametersDescr,
                             data.frame(
                               lambda = l,
                               parameter = x@regularized, 
                               min = apply(parameters[parameters$lambda == l, x@regularized],2,min), 
                               mean = apply(parameters[parameters$lambda == l, x@regularized],2,mean), 
                               max = apply(parameters[parameters$lambda == l, x@regularized],2,max))
    )
  }
  ggplot2::ggplot(parametersDescr, mapping = ggplot2::aes(x = lambda, y = mean, ymin=min, ymax=max, fill= parameter)) +
    ggplot2::geom_line()+
    ggplot2::geom_ribbon(alpha=0.2)+
    ggplot2::ggtitle("Regularized Parameters")
})
