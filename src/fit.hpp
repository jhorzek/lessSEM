#ifndef FIT_H
#define FIT_H

#include <RcppArmadillo.h>
#include "config.hpp"
// [[Rcpp :: depends ( RcppArmadillo )]]

double computeIndividualM2LL(const int nObservedVariables, const arma::colvec& rawData, const arma::colvec& impliedMeans, const arma::mat& impliedCovariance);

// Computes the -2log likelihood in a RAM model for a  group of people with identical missing structure given the sample size, number of non-missing variables (nObservedVariables),
// observedMeans within this sample, observedCov within this sample, filtered impliedMeans, and filtered impliedCovariance.
// The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05

double computeGroupM2LL(const int sampleSize, const int nObservedVariables, const arma::colvec& observedMeans, const arma::mat& observedCov,
                        const arma::colvec& impliedMeans, const arma::mat& impliedCovariance);

#endif