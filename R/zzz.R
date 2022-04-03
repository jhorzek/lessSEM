#' aCV4SEM
#'
#' aCV4SEM implements approximate cross-validation for SEM
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib aCV4SEM
#' @name aCV4SEM

Rcpp::loadModule("SEM_cpp", TRUE)
Rcpp::loadModule("SEM_cpp", TRUE)

