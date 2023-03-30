#' Class for log-likelihood of regularized SEM. Note: we define a custom logLik - 
#' Function because the generic one is using df = number of parameters which might be confusing.
#' @slot logLik log-Likelihood
#' @slot nParameters number of parameters in the model
#' @slot N number of persons in the data set
#' @export
setClass("logLikelihood",
         representation = representation(
           logLik="numeric",
           nParameters="integer", 
           N="integer"
         ))

#' show
#' 
#' @param object object of class logLikelihood
#' @return No return value, just prints estimates
#' @export
setMethod("show", "logLikelihood", function(object){
  cat(paste0(object@logLik), "(nPar = ", object@nParameters, ")")
})
