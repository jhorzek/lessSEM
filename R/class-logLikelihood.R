# Note: we define a custom logLik - Function because the generic one is 
# using df = number of parameters which might be confusing.
setClass("logLikelihood",
         representation = representation(
           logLik="numeric",
           nParameters="integer", 
           N="integer"
         ))

#' show
#' 
#' @param object object of class logLikelihood
#' @export
setMethod("show", "logLikelihood", function(object){
  cat(paste0(object@logLik), "(nPar = ", object@nParameters, ")")
})
