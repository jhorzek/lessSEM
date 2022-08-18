#ifndef SEMWITHTRANSFORMATIONS_H
#define SEMWITHTRANSFORMATIONS_H

#include "SEM.h"

// first, we have to define some types to allow for passing compiled
// functions to our SEM. These functions will be used to transform the parameters
// See: https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
typedef Rcpp::NumericVector (*transformationFunctionPtr)(
    Rcpp::NumericVector& //labeled parameter values
); // function takes our raw parameters and returns the transformed ones!
typedef Rcpp::XPtr<transformationFunctionPtr> transformationFunctionPtr_t;

class parametersWithTransformations: public parameters{
  // parametersWithTransformations is like parameters, but
  // some of its elements may be transformations from other 
  // elements -> not every parameter is a real parameter
  
public: 
  std::vector<bool> isTransformation;
  transformationFunctionPtr transformationFunction;
  
  // we use the same initialize function as that of parameters.
  // Here, we initialize all raw and transformed parameters
  
  // Now, we have to tell our SEM, which parameters are transformations
  // and how to compute those
  void asTransformation(Rcpp::StringVector parameterLabels, // labels of parameters which are
                        // transformations of other parameters
                        SEXP transformationFunctionSEXP
  );
  
  void transform();
  
};

class SEMWithTransformationsCpp: public SEMCpp{
public:
  parametersWithTransformations parameterTable;
  
  void addTransformation(Rcpp::StringVector parameterLabels, // labels of parameters which are
                         // transformations of other parameters
                         SEXP transformationFunctionSEXP
  );
  
  void computeTransformations();
  
  arma::rowvec getGradientsWithTransformation(bool raw);
  arma::mat getScoresWithTransformation(bool raw);
  arma::mat getHessianWithTransformation(Rcpp::StringVector label_,
                       arma::vec value_,
                       bool raw,
                       double eps);
  
};

RCPP_EXPOSED_CLASS(SEMWithTransformationsCpp)

#endif