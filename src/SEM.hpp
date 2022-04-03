#ifndef SEMMODULE_H
#define SEMMODULE_H

#include <RcppArmadillo.h>
#include "dataset.hpp"
#include "parameters.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class SEMCpp{
public: 
  dataset data;
  arma::mat rawData;
  Rcpp::CharacterVector manifestNames;
  
  // paramters
  parameters parameterTable;
  
  // model matrices
  arma::mat Amatrix, Smatrix, Fmatrix;
  arma::colvec mVector;
  
  // fit elements
  arma::mat impliedCovariance;
  arma::colvec impliedMeans;
  double m2LL; // minus 2 log-Likelihood
  
  // constructor
  SEMCpp(){};
  
  // setter or change elements
  void addRawData(arma::mat rawData_, Rcpp::CharacterVector manifestNames_);
  void addSubset(int N_,
                 int observed_, // number of variables without missings
                 arma::uvec notMissing_, // vector with indices of non-missing values
                 // the following elements are only relevant for N>1
                 arma::mat covariance_,
                 arma::colvec means_,
                 // raw data is required for N == 1
                 arma::colvec dataNoMissing_);
  void removeSubset(int whichSubset);
  void setMatrix(std::string whichMatrix, arma::mat values);
  void setVector(std::string whichVector, arma::colvec values);
  
  void initializeParameters(Rcpp::CharacterVector label_,
                            Rcpp::CharacterVector location_,
                            arma::uvec row_,
                            arma::uvec col_,
                            arma::vec value_,
                            arma::vec rawValue_);
  void setParameters(Rcpp::CharacterVector label_,
                     arma::vec value_,
                     bool raw);
  
  // getter
  Rcpp::DataFrame getParameters();
  
  // fit related functions
  void implied(); // compute implied means and covariance
  double fit();
  arma::colvec gradients();
  arma::mat scores();
  arma::mat hessian();
};

#endif