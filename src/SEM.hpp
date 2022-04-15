#ifndef SEMMODULE_H
#define SEMMODULE_H

#include <RcppArmadillo.h>
#include "dataset.hpp"
#include "parameters.hpp"
#include "derivativeStructure.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class SEMCpp{
public: 
  //flags
  bool wasChecked = false; // true if the model was checked
  bool wasFit = false; // true if fit was called
  int currentStatus = 0; // 0 = initialized
                     // 1 = matrices were defined
                     // 2 = initializeParameters was called
                     // 3 = addDerivativeElement was called
                     // 4 = addRawData was called
                     // 5 = addSubset was called
  
  // data
  dataset data;
  arma::mat rawData;
  arma::uvec personInSubset;
  Rcpp::StringVector manifestNames;
  
  // parameters
  parameters parameterTable;
  
  // derivative elements
  derivativeElements derivElements;
  
  // model matrices
  arma::mat Amatrix, Smatrix, Fmatrix;
  arma::colvec Mvector;
  
  // fit elements
  arma::mat impliedCovariance;
  arma::colvec impliedMeans;
  double m2LL; // minus 2 log-Likelihood
  
  // constructor
  SEMCpp(){};
  
  bool checkModel();
  
  // setter or change elements
  void addRawData(arma::mat rawData_, Rcpp::StringVector manifestNames_, arma::uvec personInSubset_);
  void addSubset(int N_,
                 int observed_, // number of variables without missings
                 arma::uvec notMissing_, // vector with indices of non-missing values
                 // the following elements are only relevant for N>1
                 arma::mat covariance_,
                 arma::colvec means_,
                 // raw data is required for N == 1
                 arma::colvec rawData_);
  void removeSubset(int whichSubset);
  void setMatrix(std::string whichMatrix, arma::mat values);
  void setVector(std::string whichVector, arma::colvec values);
  
  void initializeParameters(Rcpp::StringVector label_,
                            Rcpp::StringVector location_,
                            arma::uvec row_,
                            arma::uvec col_,
                            arma::vec value_,
                            arma::vec rawValue_);
  void setParameters(Rcpp::StringVector label_,
                     arma::vec value_,
                     bool raw);
  
  void addDerivativeElement(std::string label_, 
                            std::string location_, 
                            bool isVariance_, 
                            arma::mat positionMatrix_);
  
  // getter
  Rcpp::DataFrame getParameters();
  Rcpp::StringVector getParameterLabels();
  
  // fit related functions
  void implied(); // compute implied means and covariance
  double fit();
  arma::rowvec getGradients(bool raw);
  arma::mat getScores(bool raw);
  arma::mat getHessian(Rcpp::StringVector label_,
                       arma::vec value_,
                       bool raw,
                       double eps);
};

#endif