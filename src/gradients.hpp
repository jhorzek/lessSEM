#ifndef GRADIENTS_H
#define GRADIENTS_H

#include <RcppArmadillo.h>
#include "SEM.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::colvec computeSingleSubjectGradient_Internal(const arma::colvec& data_i, 
                                                   const SEMCpp& SEM,
                                                   const arma::mat& IminusAInverse,
                                                   const arma::uvec isObserved,
                                                   const int personInMissingGroup,
                                                   const arma::mat& impliedCovInverse,
                                                   const arma::colvec& logDetSigmas, 
                                                   const std::vector<arma::mat>& derivativesOfCovariance);

#endif