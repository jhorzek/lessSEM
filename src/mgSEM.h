#ifndef MGSEMMODULE_H
#define MGSEMMODULE_H

#include <RcppArmadillo.h>
#include "SEM.h"
#include "hessian.h"
#include "findStringInVector.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

class mgParameters{
public:
  
  arma::vec uniqueValues;
  std::vector<std::string> uniqueLabels;
  Rcpp::StringVector uniqueLabelsRcpp;
  arma::rowvec uniqueGradients;
  arma::mat uniqueHessian;
  std::vector<bool> isTransformation;
  bool hasTransformations = false;
  double gradientStepSize = 1e-6; // step size used to compute the gradients of
  // the transformations
  
  // vectors telling us where the parameters of 
  // a model are located in the parameter vector. We need two vectors of indices:
  // one which tells us where a the parameter is located in the model (e.g., position 3)
  // and one which tells us where that parameter is located in the 
  // uniqueLabels (e.g., position 1)
  std::vector<Rcpp::IntegerVector> parameterLocationInModelRcpp; // Rcpp needs integer vectors for subsetting
  std::vector<Rcpp::IntegerVector> parameterLocationInVectorRcpp; // Rcpp needs integer vectors for subsetting
  std::vector<arma::uvec> parameterLocationInModelUvec; // armadillo needs uvec for subsetting
  std::vector<arma::uvec> parameterLocationInVectorUvec; // armadillo needs uvec for subsetting
  
  std::vector<std::vector<int>> model; // saves the models where the parameters are located
  
  // in case of transformations
  transformationFunctionPtr transformationFunction;
  Rcpp::List transformationList;
  
  // we use the same initialize function as that of parameters.
  // Here, we initialize all raw and transformed parameters
  
  // Now, we have to tell our SEM, which parameters are transformations
  // and how to compute those
  void addTransformation(Rcpp::NumericVector extendedParameters,
                         std::vector<bool> isTransformation_,
                         SEXP transformationFunctionSEXP,
                         Rcpp::List transformationList_);
  
  void transform();
  
  arma::mat getTransformationGradients();
  
};



class mgSEM{
public:
  std::vector<SEMCpp> models;
  int sampleSize = 0;
  
  double m2LL;
  
  mgParameters parameters;
  arma::rowvec gradients;
  arma::mat Hessian;
  
  Rcpp::StringVector parameterLabels;
  Rcpp::NumericVector parameterValues;
  std::vector<int> locatedInModel;
  
  void addModel(Rcpp::List SEMList);
  
  void addTransformation(Rcpp::NumericVector extendedParameters,
                         std::vector<bool> isTransformation_,
                         SEXP transformationFunctionSEXP,
                         Rcpp::List transformationList_);
  
  void computeTransformations();
  
  // changing parameters
  
  void setParameters(Rcpp::StringVector label_,
                     arma::vec value_,
                     bool raw);
  
  // getter
  Rcpp::List getParameters();
  Rcpp::List getSubmodelParameters(); // full data frame including transformations
  
  Rcpp::StringVector getParameterLabels();
  
  // fit related functions
  void implied();
  bool impliedIsPD();
  double fit();
  
  arma::rowvec getGradients(bool raw);
  
  arma::mat getScores();
  
  arma::mat getHessian(Rcpp::StringVector label_,
                       arma::vec value_,
                       bool raw,
                       double eps);
  
  void setTransformationGradientStepSize(double gradientStepSize);
  
};

RCPP_EXPOSED_CLASS(mgSEM)
# endif