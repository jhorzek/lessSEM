#include <RcppArmadillo.h>
#include "implied.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

// [[Rcpp::export]]
arma::mat computeImpliedCovariance(const arma::mat& Fmatrix, const arma::mat& Amatrix, const arma::mat& Smatrix){
  
  arma::mat I = arma::eye(arma::size(Amatrix));
  arma::mat impliedCovariance = Fmatrix*arma::inv(I-Amatrix)*Smatrix*arma::trans(arma::inv(I-Amatrix))*arma::trans(Fmatrix);
  
  return(impliedCovariance);
}

// [[Rcpp::export]]
arma::mat computeImpliedMeans(const arma::mat& Fmatrix, const arma::mat& Amatrix, const arma::colvec& mVector){
  
  arma::mat I = arma::eye(arma::size(Amatrix));
  arma::mat impliedMeans = Fmatrix*arma::inv(I-Amatrix)*mVector;
  
  return(impliedMeans);
}

