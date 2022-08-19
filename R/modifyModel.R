#' modifyModel
#' 
#' Modify the model from lavaan to fit your needs
#' 
#' @param addMeans If lavaanModel has meanstructure = FALSE, addMeans = TRUE will add a mean structure. FALSE will set the means of the observed variables to the average
#' @param activeSet Option to only use a subset of the individuals in the data set. Logical vector of length N indicating which subjects should remain in the sample.
#' @param dataSet option to replace the data set in the lavaan model with a different data set. Can be useful for cross-validation
#' @param transformations allows for transformations of parameters - useful for measurement invariance tests etc.
#' @returns Object of class modifyModel
modifyModel <- function(
    addMeans = TRUE,
    activeSet = NULL,
    dataSet = NULL,
    transformations = NULL
  ){
  mod <- as.list(environment())
  class(mod) <- "modifyModel"
  return(mod)
}