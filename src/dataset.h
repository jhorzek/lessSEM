#ifndef DATAMODULE_H
#define DATAMODULE_H

#include <RcppArmadillo.h>
#include "subset.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

class dataset{
public: 
  arma::mat rawData;
  unsigned int nGroups = 0;
  std::vector<subset> dataSubsets;
  Rcpp::IntegerVector personInSubset; // tells us in which subset an individual is

  // constructor
  dataset(){};
  
  // add
  void addSubset(unsigned int N_,
                 arma::uvec persons_, // numeric; indicates which persons are in this subset
                 int observed_, // number of variables without missings
                 arma::uvec notMissing_, // vector with indices of non-missing values
                 // the following elements are only relevant for N>1
                 arma::mat covariance_,
                 arma::colvec means_,
                 // raw data is required for N == 1
                 // person in row, items in column
                 arma::mat rawData_);
  
  // remove
  void removeSubset(int whichSubset);
};

#endif