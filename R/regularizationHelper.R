# The following functions can be used to extract the 
# names of different parameters from a lavaan object.

#' loadings
#' 
#' Extract the labels of all loadings found in a lavaan model.
#' @param lavaanModel fitted lavaan model
#' @return vector with parameter labels
#' @export
#' @examples
#' # The following is adapted from ?lavaan::sem
#' library(lessSEM)
#' model <- ' 
#'   # latent variable definitions
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'   dem65 =~ y5 + a*y6 + b*y7 + c*y8
#' 
#'   # regressions
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#' 
#'   # residual correlations
#'   y1 ~~ y5
#'   y2 ~~ y4 + y6
#'   y3 ~~ y7
#'   y4 ~~ y8
#'   y6 ~~ y8
#' '
#' 
#' fit <- sem(model, data = PoliticalDemocracy)
#' 
#' loadings(fit)
loadings <- function(lavaanModel){
  
  pt <- .labelLavaanParameters(lavaanModel)
  
  return(
    pt$label[pt$op == "=~"] |>
      unique()
  )
  
}

#' regressions
#' 
#' Extract the labels of all regressions found in a lavaan model.
#' @param lavaanModel fitted lavaan model
#' @return vector with parameter labels
#' @export
#' @examples
#' # The following is adapted from ?lavaan::sem
#' library(lessSEM)
#' model <- ' 
#'   # latent variable definitions
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'   dem65 =~ y5 + a*y6 + b*y7 + c*y8
#' 
#'   # regressions
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#' 
#'   # residual correlations
#'   y1 ~~ y5
#'   y2 ~~ y4 + y6
#'   y3 ~~ y7
#'   y4 ~~ y8
#'   y6 ~~ y8
#' '
#' 
#' fit <- sem(model, data = PoliticalDemocracy)
#' 
#' regressions(fit)
regressions <- function(lavaanModel){
  pt <- .labelLavaanParameters(lavaanModel)
  
  return(
    pt$label[pt$op == "~"] |>
      unique()
  )
}

#' covariances
#' 
#' Extract the labels of all covariances found in a lavaan model.
#' @param lavaanModel fitted lavaan model
#' @return vector with parameter labels
#' @export
#' @examples
#' # The following is adapted from ?lavaan::sem
#' library(lessSEM)
#' model <- ' 
#'   # latent variable definitions
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'   dem65 =~ y5 + a*y6 + b*y7 + c*y8
#' 
#'   # regressions
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#' 
#'   # residual correlations
#'   y1 ~~ y5
#'   y2 ~~ y4 + y6
#'   y3 ~~ y7
#'   y4 ~~ y8
#'   y6 ~~ y8
#' '
#' 
#' fit <- sem(model, data = PoliticalDemocracy)
#' 
#' covariances(fit)
covariances <- function(lavaanModel){
  pt <- .labelLavaanParameters(lavaanModel)
  
  return(
    pt$label[(pt$op == "~~") & 
               (pt$lhs != pt$rhs)] |>
      unique()
  )
}

#' variances
#' 
#' Extract the labels of all variances found in a lavaan model.
#' @param lavaanModel fitted lavaan model
#' @return vector with parameter labels
#' @export
#' @examples
#' # The following is adapted from ?lavaan::sem
#' library(lessSEM)
#' model <- ' 
#'   # latent variable definitions
#'   ind60 =~ x1 + x2 + x3
#'   dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'   dem65 =~ y5 + a*y6 + b*y7 + c*y8
#' 
#'   # regressions
#'   dem60 ~ ind60
#'   dem65 ~ ind60 + dem60
#' 
#'   # residual correlations
#'   y1 ~~ y5
#'   y2 ~~ y4 + y6
#'   y3 ~~ y7
#'   y4 ~~ y8
#'   y6 ~~ y8
#' '
#' 
#' fit <- sem(model, data = PoliticalDemocracy)
#' 
#' variances(fit)
variances <- function(lavaanModel){
  pt <- .labelLavaanParameters(lavaanModel)
  
  return(
    pt$label[(pt$op == "~~") & 
               (pt$lhs == pt$rhs)] |>
      unique()
  )
}

#' .labelLavaanParameters
#' 
#' Adds labels to unlabeled parameters in the lavaan parameter table. Also 
#' removes fixed parameters.
#' @param lavaanModel fitted lavaan model
#' @return parameterTable with labeled parameters
.labelLavaanParameters <- function(lavaanModel){
  pt <- lavaan::parameterTable(lavaanModel)
  # remove fixed:
  pt <- pt[pt$free != 0,]
  # add labels when missing:
  pt$label[pt$label == ""] <- paste0(pt$lhs, pt$op, pt$rhs)[pt$label == ""]
  return(pt)
}