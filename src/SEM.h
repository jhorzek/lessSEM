#ifndef SEMMODULE_H
#define SEMMODULE_H

#include <RcppArmadillo.h>
#include "dataset.h"
#include "parameters.h"
#include "derivativeStructure.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

enum status {
  empty,
  initialized,
  changedParameters,
  filledMatrices,
  computedImplied,
  fitted
};

class SEMCpp{
public: 
  //flags
  status currentStatus = empty;
  bool wasChecked = false; // true if the model was checked
  bool wasFit = false; // true if fit was called
  bool hasTransformations = false; // true if the user defined transformations of parameters
  int functionCalls = 0;
  int gradientCalls = 0;
  
  // data
  dataset data;
  arma::mat rawData;
  arma::uvec personInSubset;
  Rcpp::StringVector manifestNames;
  int sampleSize;
  
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
  
  arma::rowvec gradients;
  arma::mat transformationGradients;
  
  // constructor
  SEMCpp(){};
  
  // fill model with elements
  void fill(Rcpp::List SEMList);
    
  bool checkModel();
  
  // setter or change elements
  void setParameters(Rcpp::StringVector label_,
                     arma::vec value_,
                     bool raw);
  
  void addTransformation(SEXP transformationFunctionSEXP,
                         Rcpp::List transformationList);
  
  void computeTransformations();
  
  // getter
  Rcpp::DataFrame getParameters();
  Rcpp::StringVector getParameterLabels();
  
  // fit related functions
  void implied(); // compute implied means and covariance
  bool impliedIsPD(); // check if model implied covariance matrix is positive definite
  double fit();
  arma::rowvec getGradients(bool raw);
  arma::mat getScores(bool raw);
  arma::mat getHessian(Rcpp::StringVector label_,
                       arma::vec value_,
                       bool raw,
                       double eps);
  
};

RCPP_EXPOSED_CLASS(SEMCpp)
  
#endif