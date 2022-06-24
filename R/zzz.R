#' linr
#'
#' linr implements approximate cross-validation for SEM
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib linr
#' @name linr

Rcpp::loadModule("SEM_cpp", TRUE)
Rcpp::loadModule("istaEnet_cpp", TRUE)
Rcpp::loadModule("istaEnetGeneralPurpose_cpp", TRUE)