#include <RcppArmadillo.h>
#include "bfgs.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

// [[Rcpp::export]]
arma::mat BFGS_C(const arma::rowvec& parameters_kMinus1, 
                 const arma::rowvec& gradients_kMinus1, 
                 const arma::mat& Hessian_kMinus1, 
                 const arma::rowvec& parameters_k, 
                 const arma::rowvec& gradients_k, 
                 const bool cautious, 
                 const double hessianEps){
  
  return(lessSEM::BFGS(
    parameters_kMinus1, 
    gradients_kMinus1, 
    Hessian_kMinus1, 
    parameters_k, 
    gradients_k, 
    cautious, 
    hessianEps
  ));
}