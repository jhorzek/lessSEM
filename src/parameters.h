#ifndef PARAMETERMODULE_H
#define PARAMETERMODULE_H

#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp :: depends ( RcppArmadillo )]]
namespace parameterModule{
struct parameterElements{
  // every parameter value can only have one value
  // but may appear in multiple rows and columns of a matrix
  double rawValue; // raw (untransformed) parameter value
  double value; // parameter value
  bool changed; // to indicate if the parameter changed since the last fit
  bool isVariance;
  std::string location; // name of the matrix where the parameter is located
  bool isTransformation;
  
  std::vector<int> row;
  std::vector<int> col;
  //arma::uvec row = {0}; // row in the matrix
  //arma::uvec col = {0}; // col in the matrix
};
}

// In some SEMs, we also allow for parameters to be transformations of one another.
// To this end, we define some types to allow for passing compiled
// functions to our SEM. These functions will be used to transform the parameters
// See: https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
typedef Rcpp::NumericVector (*transformationFunctionPtr)(
    Rcpp::NumericVector&, //labeled parameter values
    Rcpp::List& //additional arguments
); // function takes our raw parameters and returns the transformed ones!
typedef Rcpp::XPtr<transformationFunctionPtr> transformationFunctionPtr_t;


class parameters{
public: 
  
  std::map< std::string, parameterModule::parameterElements> parameterMap;
  
  // for export to R
  Rcpp::CharacterVector uniqueParameterLabels;
  Rcpp::NumericVector uniqueParameterValues;
  Rcpp::NumericVector uniqueRawParameterValues;
  Rcpp::CharacterVector uniqueParameterLocations;
  std::vector<bool> uniqueIsTransformation;
  
  bool hasTransformations = false;
  double gradientStepSize = 1e-6; // step size used to compute the gradients of
  // the transformations
  
  // to decide which elements should be recomputed
  bool AChanged = true;
  bool SChanged = true;
  bool mChanged = true;
  
  int nModelParameters, nRealParameters;
  // nModelParameters: number of parameters in the SEM, some of which
  // can be functions of other parameters defined in the transformation
  // nRealParameters: number of actual parameters, excluding 
  // those that are transformations
  
  // constructor
  parameters(){}
  
  //
  void initialize(Rcpp::StringVector label_,
                  Rcpp::StringVector location_,
                  arma::uvec row_,
                  arma::uvec col_,
                  arma::vec value_,
                  arma::vec rawValue_,
                  std::vector<bool> isTransformation_);
  
  // getter
  Rcpp::DataFrame getParameters();
  
  // setter
  void setParameters(Rcpp::StringVector label_,
                     arma::vec value_,
                     bool raw
  );
  
  // in case of transformations
  transformationFunctionPtr transformationFunction;
  Rcpp::List transformationList;
  
  // we use the same initialize function as that of parameters.
  // Here, we initialize all raw and transformed parameters
  
  // Now, we have to tell our SEM, which parameters are transformations
  // and how to compute those
  void addTransformation(SEXP transformationFunctionSEXP,
                         Rcpp::List transformationList_);
  
  void transform();
  
  arma::mat getTransformationGradients();
  
};

#endif