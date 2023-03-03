#ifndef WLS_h
#define WLS_h
#include <RcppArmadillo.h>

double WLS(const arma::mat& weightsInverse,
           const arma::mat& observedCov,
           const arma::colvec& observedMeans,
           const arma::colvec& impliedMeans, 
           const arma::mat& impliedCovariance);

double WLSDerivative(const arma::mat& weightsInverse,
                     const arma::mat& observedCov,
                     const arma::colvec& observedMeans,
                     const arma::colvec& impliedMeans, 
                     const arma::mat& impliedCovariance,
                     const arma::colvec& impliedMeansDerivative,
                     const arma::mat& impliedCovarianceDerivative);
#endif