setClass(Class = "aCV4regularizedSEMWithCustomPenalty",
         representation = representation(
           parameters="data.frame",
           tuningParameters="data.frame",
           cvfits = "data.frame",
           parameterLabels = "character",
           cvfitsDetails="data.frame", 
           subsets = "matrix"
         )
)

#' show
#' 
#' @export
setMethod("show", "aCV4regularizedSEMWithCustomPenalty", function (object) {
  bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
  cat(paste0("Best parameter estimates observed for: ", paste0(names(object@tuningParameters), " = ", object@tuningParameters[bestFit,]), collapse = ", "), "\n\n")
  
  parameters <- unlist(object@parameters[bestFit,object@parameterLabels])
  print(parameters)
})

#' summary
#' 
#' @export
setMethod("summary", "aCV4regularizedSEMWithCustomPenalty", function (object) {
  modelName <-deparse(substitute(object)) # get the name of the object
  cat(paste0("#### Approximate Cross Validation Results ####\n\n"))
  
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
#' Returns the parameter estimates of an aCV4regularizedSEMWithCustomPenalty
#' 
#' @param object object of class regularizedSEM
#' 
#' @export
setMethod("coef", "aCV4regularizedSEMWithCustomPenalty", function (object) {
  
  bestFit <- which(object@cvfits$cvfit == min(object@cvfits$cvfit))
  parameters <- unlist(object@parameters[bestFit,object@parameterLabels])
  return(parameters)
  
})

#' plot
#' 
#' plots the cross-validation fits
#' 
#' @param x object of class aCV4regularizedSEMWithCustomPenalty
#' @export
setMethod("plot", "aCV4regularizedSEMWithCustomPenalty", function (x) {
  
  parameters <- x@parameters
  fits <- x@cvfits

  plot(x = 1:nrow(fits), 
       y = fits$cvfit, 
       type = "l",
       ylab = "cv-fit", 
       xlab = "index",
       main = "index refers to the row in object@tuningParameters[row,].")
  points(x = 1:nrow(fits), 
         y = fits$cvfit)
})