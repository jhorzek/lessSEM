#' S4 method to compute fit indices (e.g., AIC, BIC, ...)
#'
#' @param object a model fitted with lessSEM
#' @return returns a data.frame with fit indices
#' @export
setGeneric("fitIndices", function(object) {
  standardGeneric("fitIndices")
})

#' S4 method to exract the estimates of an object
#'
#' @param object a model fitted with lessSEM
#' @param criterion fitIndice used to select the parameters
#' @param transformations boolean: Should transformations be returned?
#' @return returns a matrix with estimates
#' @export
setGeneric("estimates", function(object, criterion = NULL, transformations = FALSE) {
  standardGeneric("estimates")
})
