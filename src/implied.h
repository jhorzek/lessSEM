#ifndef IMPLIED_H
#define IMPLIED_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat computeImpliedCovarianceFull(const arma::mat& Amatrix,
                                       const arma::mat& Smatrix,
                                       const arma::mat& IminusAInverse);

arma::mat computeImpliedCovariance(const arma::mat& impliedCovarianceFull,
                                   const arma::mat& Smatrix);

arma::colvec computeImpliedMeansFull(const arma::mat& Amatrix, 
                                     const arma::colvec& Mvector,
                                     const arma::mat& IminusAInverse);

arma::colvec computeImpliedMeans(const arma::mat& Fmatrix,
                                 const arma::mat& impliedMeansFull);

#endif