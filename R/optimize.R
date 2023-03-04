#' bfgs
#' 
#' This function allows for optimizing models built in lavaan using the BFGS optimizer
#' implemented in lessSEM. Its elements can be accessed
#' with the "@" operator (see examples). The main purpose is to make transformations
#' of lavaan models more accessible.
#' 
#' @param lavaanModel model of class lavaan 
#' @param modifyModel used to modify the lavaanModel. See ?modifyModel.
#' @param control used to control the optimizer. See ?controlBFGS for more details.
#' @returns Model of class regularizedSEM
#' @examples 
#' library(lessSEM)
#' 
#' # Identical to regsem, lessSEM builds on the lavaan
#' # package for model specification. The first step
#' # therefore is to implement the model in lavaan.
#' 
#' dataset <- simulateExampleData()
#' 
#' lavaanSyntax <- "
#' f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 +
#'      l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 +
#'      l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
#' f ~~ 1*f
#' "
#' 
#' lavaanModel <- lavaan::sem(lavaanSyntax,
#'                            data = dataset,
#'                            meanstructure = TRUE,
#'                            std.lv = TRUE)
#' 
#' 
#' lsem <- bfgs(
#'   # pass the fitted lavaan model
#'   lavaanModel = lavaanModel)
#' 
#' # the coefficients can be accessed with:
#' coef(lsem)
#' 
#' # elements of lsem can be accessed with the @ operator:
#' lsem@parameters
#' @export
bfgs <- function(lavaanModel,
                 modifyModel = lessSEM::modifyModel(),
                 control = lessSEM::controlBFGS()){
  
  regularized <- NULL
  
  result <- .regularizeSmoothSEMInternal(
    lavaanModel = lavaanModel,
    penalty = "ridge",
    weights = regularized,
    tuningParameters = data.frame(lambda = 0,
                                  alpha = 0),
    epsilon = 0, # ridge is already smooth
    tau = 0,
    modifyModel = modifyModel,
    control = control
  )
  return(result)
  
}