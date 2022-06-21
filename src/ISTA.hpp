#ifndef SEMISTA_h
#define SEMISTA_h
#include "fista.h"
#include "SEM.hpp"

// [[Rcpp :: depends ( fista )]]

class SEMIsta: public ista::model{
public:
  
  SEMCpp SEM;
  
  SEMIsta(SEMCpp SEM_){
    SEM = SEM_;
  }
  
  void setParameters(arma::rowvec parameterValues) override {
    SEM.setParameters(SEM.getParameterLabels(), // labels
                      arma::trans(parameterValues), // values
                      true); // raw
  }
  
  double fit() override{
    try{
    SEM.fit();
    }catch(...){
      Rcpp::Rcout << "Fit error" << std::endl;
      Rcpp::Rcout << SEM.m2LL << std::endl;
      return(arma::datum::nan);
    }
    if(!SEM.impliedCovariance.is_sympd()){
      Rcpp::Rcout << "not positive definite" << std::endl;
      return(arma::datum::nan);
    }
    return(SEM.m2LL);
  }
  
  arma::rowvec gradients() override{
    return(SEM.getGradients(true)); 
  }
  
};

#endif