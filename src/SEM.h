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

enum estimator{
  fiml,
  wls
};

struct WLS{
public:
  arma::mat weightsInverse;
  bool meanstructure;
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
  estimator estim;
  
  WLS WLSElements;
  
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
  arma::colvec impliedMeans;
  arma::colvec impliedMeansFull;
  arma::mat impliedCovariance;
  arma::mat impliedCovarianceFull;
  arma::mat IminusAInverse; // this element is used repeatedly, so we
  // define it once and reuse it.
  double objectiveValue; // minus 2 log-Likelihood
  
  // in case of missing data, we have to compute
  // means and covariances for each subset. These are going
  // to be saved in the following vectors:
  std::vector<arma::colvec> subsetImpliedMeans;
  std::vector<arma::mat> subsetImpliedCovariance;
  std::vector<arma::mat> subsetImpliedCovarianceInverse;
  
  arma::rowvec gradients;
  arma::mat transformationGradients;
  
  // We have to compute
  // derivatives for the means and covariances with respect to the parameters. 
  // These are going to be saved in the following vectors:
  bool detivativesInitialized = false;
  std::vector<arma::mat> impliedCovarianceDerivatives;
  std::vector<arma::colvec> impliedMeansDerivatives;
  
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
  
  std::string getEstimator();
  
  arma::rowvec getGradients(bool raw);
  arma::mat getScores(bool raw);
  arma::mat getHessian(Rcpp::StringVector label_,
                       arma::vec value_,
                       bool raw,
                       double eps);
  
  void setTransformationGradientStepSize(double gradientStepSize);
  
};

RCPP_EXPOSED_CLASS(SEMCpp)
  
#endif