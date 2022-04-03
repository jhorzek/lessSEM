#ifndef IMPLIED_H
#define IMPLIED_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat computeImpliedCovariance(const arma::mat& Fmatrix, 
                                   const arma::mat& Amatrix,
                                   const arma::mat& Smatrix);

arma::mat computeImpliedMeans(const arma::mat& Fmatrix,
                              const arma::mat& Amatrix, 
                              const arma::colvec& mVector);

#endif