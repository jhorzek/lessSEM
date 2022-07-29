#ifndef SEMMODULE_H
#define SEMMODULE_H

#include <RcppArmadillo.h>
#include "dataset.hpp"
#include "parameters.hpp"
#include "derivativeStructure.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

enum status {
  addedMatrices,
  addedParameters,
  addedDerivatives,
  addedRawData,
  addedSubsets,
  changedParameters,
  filledMatrices,
  computedImplied,
  fitted
};

//' @name SEMCpp
//' 
//' @title internal SEM representation
//' 
//' @field new Creates a new SEMCpp.
//' @field A Matrix with directed effects
//' @field S Matrix with undirected paths
//' @field F Filter matrix to separate latent and observed variables
//' @field m Vector with means of observed and latent variables
//' @field impliedCovariance implied covariance matrix
//' @field impliedMeans implied means vector
//' @field m2LL -2 log-likelihood
//' @field rawData raw data set
//' @field manifestNames names of manifest variables
//' @field wasFit names of manifest variables
//' @field setMatrix Fills the elements of a model matrix. Expects a char 
//' (A, S, or F), and a matrix with values
//' @field setVector Fills the elements of a model vector. 
//' Expects a char (m), and a vector with values
//' @field initializeParameters Initializes the parameters of the model. 
//' Expects StringVector with labels, StringVector with location, 
//' uvec with rows, uvec with columns, vec with values, and vec with rawValues
//' @field setParameters Changes the parameters of a model. Expects StringVector 
//' with labels, vec with values, and boolean indicating if the values are raw (TRUE) 
//' or transformed (FALSE).
//' @field addDerivativeElement Add an element to the internal derivative structure. 
//' string label, string with location, and boolean indicating if it is a variance, matrix with positions.
//' @field addRawData Adds a raw data set. Expects matrix and vector with labels of manifest variables.
//' @field addSubset Adds a subset to the data. Expects the sample size N of the subset, 
//' the number of observed variables without missings, uvec with indices of notMissing values, 
//' matrix with covariance, colvec with means, raw data without missings.
//' @field implied Computes implied means and covariance matrix
//' @field fit Fits the model. Returns -2 log likelihood
//' @field getParameters Returns a data frame with model parameters.
//' @field getParameterLabels Returns a vector with unique parameter labels as used internally.
//' @field getGradients Returns a matrix with scores.
//' @field getScores Returns a matrix with scores.
//' @field getHessian Returns the hessian of the model. Expects the labels of the 
//' parameters and the values of the parameters as well as a boolean indicating if 
//' these are raw. Finally, a double (eps) controls the precision of the approximation.

class SEMCpp{
public: 
  //flags
  bool wasChecked = false; // true if the model was checked
  bool wasFit = false; // true if fit was called
  
  status currentStatus;
  
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
                 arma::uvec persons_, // indices for which persons are in this subset
                 int observed_, // number of variables without missings
                 arma::uvec notMissing_, // vector with indices of non-missing values
                 // the following elements are only relevant for N>1
                 arma::mat covariance_,
                 arma::colvec means_,
                 // raw data is required for N == 1
                 arma::mat rawData_);
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

RCPP_EXPOSED_CLASS(SEMCpp)
  
#endif