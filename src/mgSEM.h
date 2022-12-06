#ifndef MGSEMMODULE_H
#define MGSEMMODULE_H

#include <RcppArmadillo.h>
#include "SEM.h"
#include "hessian.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

int findStringInVector(std::string what, std::vector<std::string> where, bool throwError){
  
  for (int i = 0; i < where.size(); i++){
    if(where.at(i).compare(what) == 0){
      return(i);
    }
  }
  if(throwError){
    Rcpp::stop("Could not find parameter.");
  }else{
    return -1;
  }
}

int findStringInVector(std::string what, Rcpp::StringVector where, bool throwError){
  std::string currentString;
  for (int i = 0; i < where.size(); i++){
    currentString = Rcpp::as<std::string>(where.at(i));
    if(currentString.compare(what) == 0){
      return(i);
    }
  }
  if(throwError){
    Rcpp::stop("Could not find parameter.");
  }else{
    return -1;
  }
}

class mgParameters{
public:
  
  arma::vec uniqueValues;
  std::vector<std::string> uniqueLabels;
  Rcpp::StringVector uniqueLabelsRcpp;
  arma::colvec uniqueGradients;
  arma::mat uniqueHessian;
  std::vector<bool> isTransformation;
  bool hasTransformations = false;
  
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
  
  // vectors with values and labels for each model
  std::vector<Rcpp::NumericVector> modelParameterValues;
  std::vector<Rcpp::StringVector> modelParameterLabels;
  
  // in case of transformations
  transformationFunctionPtr transformationFunction;
  Rcpp::List transformationList;
  
  // we use the same initialize function as that of parameters.
  // Here, we initialize all raw and transformed parameters
  
  // Now, we have to tell our SEM, which parameters are transformations
  // and how to compute those
  void addTransformation(
      std::vector<bool> isTransformation_,
      SEXP transformationFunctionSEXP,
      Rcpp::List transformationList_);
  
  void transform();
  
  arma::mat getTransformationGradients();
  
};



class mgSEM{
public:
  std::vector<SEMCpp> models;
  
  mgParameters parameters;
  arma::rowvec gradients;
  arma::mat Hessian;
  
  Rcpp::StringVector parameterLabels;
  Rcpp::NumericVector parameterValues;
  std::vector<int> locatedInModel;
  
  void addModel(Rcpp::List SEMList);
  
  void addTransformation(std::vector<bool> isTransformation_,
                         SEXP transformationFunctionSEXP,
                         Rcpp::List transformationList);
  
  void computeTransformations();
  
  // changing parameters
  
  void setParameters(Rcpp::StringVector label_,
                     arma::vec value_,
                     bool raw);
  
  // getter
  Rcpp::NumericVector getParameters();
  
  Rcpp::StringVector getParameterLabels();
  
  // fit related functions
  void implied();
  
  double fit();
  
  arma::rowvec getGradients(bool t);
  
  arma::mat getScores();
  
  arma::mat getHessian(Rcpp::StringVector label_,
                       arma::vec value_,
                       double eps);
  
};

# endif