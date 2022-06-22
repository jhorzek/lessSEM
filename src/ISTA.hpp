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
  
  double fit(Rcpp::NumericVector parameterValues) override{
    
    try{
      // change parameter values
      Rcpp::StringVector parameterLabels = parameterValues.names();
      SEM.setParameters(parameterLabels, // labels
                        parameterValues, // values
                        true); // raw
      SEM.fit();
      
    }catch(...){
      return(arma::datum::nan);
    }
    if(!SEM.impliedCovariance.is_sympd()){
      return(arma::datum::nan);
    }
    return(SEM.m2LL);
  }
  
  Rcpp::NumericVector gradients(Rcpp::NumericVector parameterValues) override{
    
    Rcpp::NumericVector gradients(parameterValues.length());
    
    try{
      // change parameter values
      Rcpp::StringVector parameterLabels = parameterValues.names();
      SEM.setParameters(parameterLabels, // labels
                        parameterValues, // values
                        true); // raw
      SEM.fit();
      arma::rowvec gradients_arma = SEM.getGradients(true);
      gradients = Rcpp::NumericVector(gradients_arma.begin(), 
                                                          gradients_arma.end());
    }catch(...){
      gradients.fill(NA_REAL);
      return(gradients);
    }
    if(!SEM.impliedCovariance.is_sympd()){
      gradients.fill(NA_REAL);
      return(gradients);
    }
    
    return(gradients); 

  }
  
};

#endif