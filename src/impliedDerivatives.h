#ifndef IMPLIED_DERIVATIVES_H
#define IMPLIED_DERIVATIVES_H

#include <RcppArmadillo.h>

// the following struct will hold elements that need to be computed repeatedly:
struct derivPrecompute{
  arma::mat FIminusAInverse;
  arma::mat tFIminusAInverse;
  arma::mat FimpliedCovarianceFull;
  arma::mat impliedCovarianceFulltF;
};

arma::mat impliedCovarianceDerivative(const double parameterValue, // only required for variances due to the internal transformation used by lessSEM.
                                      const std::string& location, 
                                      const bool isVariance, // only required for variances due to the internal transformation used by lessSEM.
                                      const bool raw, // only required for variances due to the internal transformation used by lessSEM.
                                      const arma::mat& impliedCovariance, 
                                      const derivPrecompute& precomputedElements,
                                      const arma::mat& Fmatrix,
                                      const arma::mat& IminusAInverse,
                                      const arma::mat& derivativeElement);

arma::colvec impliedMeansDerivative(const std::string& location, 
                                    const arma::mat& impliedMeans, 
                                    const arma::mat& impliedMeansFull,
                                    const arma::mat& Fmatrix,
                                    const derivPrecompute& precomputedElements,
                                    const arma::mat& IminusAInverse,
                                    const arma::mat& derivativeElement);

#endif