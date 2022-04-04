#ifndef SCORES_H
#define SCORES_H

#include <RcppArmadillo.h>
#include "SEM.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat scores(const SEMCpp& SEM, bool raw);

std::vector<arma::mat> computeimpliedCovarianceDerivatives(const SEMCpp& SEM,
                                                           const arma::mat& IminusAInverse,
                                                           bool raw);
arma::mat computeLogDetSigmas(const SEMCpp& SEM,
                              const arma::mat& IminusAInverse,
                              const std::vector<arma::mat>& derivativesOfCovariance,
                              const std::vector<arma::mat>& impliedCovInverses,
                              bool raw);

#endif