#' modifyModel
#' 
#' Modify the model from lavaan to fit your needs
#' 
#' @param addMeans If lavaanModel has meanstructure = FALSE, addMeans = TRUE will add a mean structure. FALSE will set the means of the observed variables to their observed means.
#' @param activeSet Option to only use a subset of the individuals in the data set. Logical vector of length N indicating which subjects should remain in the sample.
#' @param dataSet option to replace the data set in the lavaan model with a different data set. Can be useful for cross-validation
#' @param transformations allows for transformations of parameters - useful for measurement invariance tests etc.
#' @param transformationList optional list used within the transformations. NOTE: This must be used as an Rcpp::List.
#' @param transformationGradientStepSize step size used to compute the gradients of the
#' transformations
#' @examples
#' modification <- modifyModel(addMeans = TRUE) # adds intercepts to a lavaan object
#' # that was fitted without explicit intercepts
#' @returns Object of class modifyModel
#' @export
modifyModel <- function(
    addMeans = FALSE,
    activeSet = NULL,
    dataSet = NULL,
    transformations = NULL,
    transformationList = list(),
    transformationGradientStepSize = 1e-6
  ){
  mod <- as.list(environment())
  class(mod) <- "modifyModel"
  return(mod)
}