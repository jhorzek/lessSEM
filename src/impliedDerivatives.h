#ifndef IMPLIED_DERIVATIVES_H
#define IMPLIED_DERIVATIVES_H

#include <RcppArmadillo.h>

arma::mat impliedCovarianceDerivative(const std::string& location, 
                                      const arma::mat& impliedCovariance, 
                                      const arma::mat& impliedCovarianceFull,
                                      const arma::mat& Fmatrix,
                                      const arma::mat& IminusAInverse,
                                      const arma::mat& derivativeElement);

arma::mat impliedMeansDerivative(const std::string& location, 
                                 const arma::mat& impliedMeans, 
                                 const arma::mat& impliedMeansFull,
                                 const arma::mat& Fmatrix,
                                 const arma::mat& IminusAInverse,
                                 const arma::mat& derivativeElement);

#endif