#include <RcppArmadillo.h>
#include "subset.h"
#include "dataset.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

// add

void dataset::addSubset(unsigned int N_,
                        arma::uvec persons_, // indices for which persons are in this subset
                        int observed_, // number of variables without missings
                        arma::uvec notMissing_, // vector with indices of non-missing values
                        // the following elements are only relevant for N>1
                        arma::mat covariance_,
                        arma::colvec means_,
                        // raw data is required for N == 1
                        // person in row, items in column
                        arma::mat rawData_){
  if(rawData_.n_rows != N_){
    Rcpp::stop("The number of rows of rawData_ does not match N_ in addSubset");
  }
  
  subset newSubset;
  newSubset.N = N_;
  newSubset.persons = persons_;
  newSubset.observed = observed_;
  newSubset.notMissing = notMissing_;
  newSubset.covariance = covariance_;
  newSubset.means = means_;
  newSubset.rawData = rawData_;
  newSubset.objectiveValue = 0.0;
  
  dataSubsets.push_back(newSubset);

  nGroups += 1;

  return;
}

// remove
void dataset::removeSubset(int whichSubset){
  
  dataSubsets.erase(dataSubsets.begin() + whichSubset);
  
  return;
}
