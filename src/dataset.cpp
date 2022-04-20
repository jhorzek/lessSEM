#include <RcppArmadillo.h>
#include "config.hpp"
#include "subset.hpp"
#include "dataset.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

// add

void dataset::addSubset(int N_,
                        int observed_, // number of variables without missings
                        arma::uvec notMissing_, // vector with indices of non-missing values
                        // the following elements are only relevant for N>1
                        arma::mat covariance_,
                        arma::colvec means_,
                        // raw data is required for N == 1
                        arma::colvec rawData_){
  subset newSubset;
  newSubset.N = N_;
  newSubset.observed = observed_;
  newSubset.notMissing = notMissing_;
  newSubset.covariance = covariance_;
  newSubset.means = means_;
  newSubset.rawData = rawData_;
  newSubset.m2LL = 0.0;
  
  dataSubsets.push_back(newSubset);

  nGroups += 1;

  return;
}

// remove
void dataset::removeSubset(int whichSubset){
  
  dataSubsets.erase(dataSubsets.begin() + whichSubset);
  
  return;
}
