#ifndef WLS_h
#define WLS_h
#include <RcppArmadillo.h>

// Weighted Least Squares when no means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @return double indicating fit 
double WLS(const arma::mat& weightsInverse,
           const arma::mat& observedCov,
           const arma::mat& impliedCovariance);

// Weighted Least Squares when means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @parma observedMeans observed means
// @param impliedMeans model implied means
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @return double indicating fit 
double WLS(const arma::mat& weightsInverse,
           const arma::colvec& observedMeans,
           const arma::colvec& impliedMeans, 
           const arma::mat& impliedCovariance,
           const arma::mat& observedCov);

// WLSDerivative Least Squares when no means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @param impliedCovarianceDerivative derivative of covariance matrix wrt parameter
// @return double indicating derivative of functions wrt individual parameter 
double WLSDerivative(const arma::mat& weightsInverse,
                     const arma::mat& observedCov,
                     const arma::mat& impliedCovariance,
                     const arma::mat& impliedCovarianceDerivative);

// WLSDerivative Least Squares when means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @parma observedMeans observed means
// @param impliedMeans model implied means
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @param impliedCovarianceDerivative derivative of covariance matrix wrt parameter
// @return double indicating derivative of functions wrt individual parameter 
double WLSDerivative(const arma::mat& weightsInverse,
                     const arma::mat& observedCov,
                     const arma::colvec& observedMeans,
                     const arma::colvec& impliedMeans, 
                     const arma::mat& impliedCovariance,
                     const arma::colvec& impliedMeansDerivative,
                     const arma::mat& impliedCovarianceDerivative);
#endif