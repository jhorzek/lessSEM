#ifndef SUBSETMODULE_H
#define SUBSETMODULE_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

class subset{
public: 
  int N; // number of persons in this subset
  arma::uvec persons;
  int observed; // number of variables without missings
  arma::uvec notMissing; // vector with indices of non-missing values
  // the following elements are only relevant for N>1
  arma::mat covariance;
  arma::colvec means;
  // data (with missings) is required for N == 1
  arma::mat rawData;
  double m2LL; // minus 2 log-Likelihood for this subset
  
  // constructor
  subset(){};
};

#endif