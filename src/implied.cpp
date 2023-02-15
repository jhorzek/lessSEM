#include <RcppArmadillo.h>
#include "implied.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeImpliedCovariance
//' 
//' computes the implied covariance matrix of a SEM using RAM notation
//' @param Fmatrix filter matrix
//' @param Amatrix matrix with directed effects
//' @param Smatrix matrix with undirected effects
//' @returns matrix with implied covariances
// [[Rcpp::export]]
arma::mat computeImpliedCovariance(const arma::mat& Fmatrix, const arma::mat& Amatrix, const arma::mat& Smatrix){
  
  arma::mat IminusAInverse = arma::inv(arma::eye(arma::size(Amatrix)) - Amatrix);
  arma::mat impliedCovariance = Fmatrix*IminusAInverse*Smatrix*arma::trans(IminusAInverse)*arma::trans(Fmatrix);
  
  return(impliedCovariance);
}

//' computeImpliedMeans
//' 
//' computes the implied means vector of a SEM using RAM notation
//' @param Fmatrix filter matrix
//' @param Amatrix matrix with directed effects
//' @param Mvector vector with means
//' @returns matrix with implied means
// [[Rcpp::export]]
arma::mat computeImpliedMeans(const arma::mat& Fmatrix, const arma::mat& Amatrix, const arma::colvec& Mvector){
  
  arma::mat I = arma::eye(arma::size(Amatrix));
  arma::mat impliedMeans = Fmatrix*arma::solve(I-Amatrix, Mvector);
  
  return(impliedMeans);
}

