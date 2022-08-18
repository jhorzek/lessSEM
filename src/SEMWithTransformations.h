#ifndef SEMWITHTRANSFORMATIONS_H
#define SEMWITHTRANSFORMATIONS_H
// [[Rcpp :: depends ( RcppArmadillo )]]

#include "SEM.h"

class SEMWithTransformationsCpp: public SEMCpp{
public:
double mynewvariable = 1.0;
  // void addTransformation(SEXP transformationFunctionSEXP);
  // 
  // void computeTransformations();
  
  // arma::rowvec getGradientsWithTransformation(bool raw);
  // arma::mat getScoresWithTransformation(bool raw);
  // arma::mat getHessianWithTransformation(Rcpp::StringVector label_,
  //                      arma::vec value_,
  //                      bool raw,
  //                      double eps);
  
};

RCPP_EXPOSED_CLASS(SEMWithTransformationsCpp)

#endif