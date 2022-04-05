#ifndef GRADIENTS_H
#define GRADIENTS_H

#include <RcppArmadillo.h>
#include "SEM.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::rowvec gradientsByGroup(const SEMCpp& SEM, bool raw);
  
arma::colvec computeSingleSubjectGradient_Internal(const arma::colvec& data_i, 
                                                   const SEMCpp& SEM,
                                                   const arma::mat& IminusAInverse,
                                                   const arma::uvec isObserved,
                                                   const int personInSubset,
                                                   const arma::mat& impliedCovInverse,
                                                   const arma::colvec& logDetSigmas, 
                                                   const std::vector<arma::mat>& derivativesOfCovariance);

arma::colvec computeSubgroupGradient_Internal(const int N,
                                              const arma::colvec& means, 
                                              const arma::mat& covariance,
                                              const SEMCpp& SEM,
                                              const arma::mat& IminusAInverse,
                                              const arma::uvec isObserved,
                                              const int groupInSubset,
                                              const arma::mat& impliedCovInverse,
                                              const arma::colvec& logDetSigmas, 
                                              const std::vector<arma::mat>& derivativesOfCovariance);

#endif