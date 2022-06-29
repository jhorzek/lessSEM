#ifndef SEMISTA_h
#define SEMISTA_h
#include "SEM.hpp"
#include "lessSEM.hpp"

class SEMFitFramework: public lessSEM::model{
public:
  
  SEMCpp& SEM;
  
  SEMFitFramework(SEMCpp& SEM_): SEM(SEM_){}
    
  double fit(arma::rowvec parameterValues,
             Rcpp::StringVector parameterLabels) override{
    
    try{
      // change parameter values
      SEM.setParameters(parameterLabels, // labels
                        parameterValues.t(), // values
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
  
  arma::rowvec gradients(arma::rowvec parameterValues,
                         Rcpp::StringVector parameterLabels) override{
    
    arma::rowvec gradients(parameterValues.n_elem);
    
    try{
      // change parameter values
      SEM.setParameters(parameterLabels, // labels
                        parameterValues.t(), // values
                        true); // raw
      SEM.fit();
      gradients = SEM.getGradients(true);
      
    }catch(...){
      gradients.fill(arma::datum::nan);
      return(gradients);
    }
    if(!SEM.impliedCovariance.is_sympd()){
      gradients.fill(arma::datum::nan);
      return(gradients);
    }
    
    return(gradients); 

  }
  
};

#endif