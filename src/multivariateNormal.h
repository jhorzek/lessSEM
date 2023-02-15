#ifndef MVNORM_h
#define MVNORM_h
#include <RcppArmadillo.h>

double m2LLMultiVariateNormal(const arma::colvec& observed,
                              const arma::colvec& impliedMeans, 
                              const arma::mat& impliedCovariance,
                              const arma::mat& impliedCovarianceInverse,
                              const double logDetImpliedCovariance);

double m2LLMultiVariateNormalDerivative( const std::string& location,
                                         const arma::colvec& observed,
                                         const arma::mat& impliedMeans,
                                         const arma::colvec& impliedMeansDerivative,
                                         const arma::mat& impliedCovariance, 
                                         const arma::mat& impliedCovarianceInverse,
                                         const arma::mat& impliedCovarianceDerivative);

double m2LLGroupMultiVariateNormal(double N,
                                   const arma::colvec& observedMeans, 
                                   const arma::colvec& impliedMeans, 
                                   const arma::mat& observedCov,
                                   const arma::mat& impliedCovariance,
                                   const arma::mat& impliedCovarianceInverse,
                                   const double logDetImpliedCovariance);

double m2LLGroupMultiVariateNormalDerivative(
    const std::string& location,
    double N,
    const arma::colvec& observedMeans, 
    const arma::colvec& impliedMeans, 
    const arma::colvec& impliedMeansDerivative,
    const arma::mat& observedCov,
    const arma::mat& impliedCovariance,
    const arma::mat& impliedCovarianceInverse,
    const arma::mat& impliedCovarianceDerivative
);

#endif