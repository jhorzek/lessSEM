#' Class for the coefficients estimated by lessSEM.
#' @slot tuningParameters tuning parameters
#' @slot estimates parameter estimates
setClass("lessSEMCoef",
         representation = representation(
           tuningParameters = "data.frame",
           estimates = "matrix"
         ))

#' show
#' 
#' @param object object of class lessSEMCoef
setMethod("show", "lessSEMCoef", function (object) {
  cat("Parameter estimates:\n")
  print(cbind(object@tuningParameters, as.data.frame(object@estimates)))
})