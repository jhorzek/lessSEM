#ifndef GPFITFRAMEWORK_h
#define GPFITFRAMEWORK_h
#include "SEM.hpp"
#include "linr.hpp"

class generalPurposeFitFramework: public linr::model{
public:
  // functions
  Rcpp::Function fitFunction;
  Rcpp::Function gradientFunction;
  Rcpp::List userSuppliedElements;
  
  generalPurposeFitFramework(
    Rcpp::Function fitFunction_,
    Rcpp::Function gradientFunction_,
    Rcpp::List userSuppliedElements_
  ): fitFunction(fitFunction_),
  gradientFunction(gradientFunction_),
  userSuppliedElements(userSuppliedElements_){}
  
  double fit(arma::rowvec parameterValues,
             Rcpp::StringVector parameterLabels) override{
    double fitValue;
    try{
      // compute fit
      fitValue = Rcpp::as<double>(fitFunction(parameterValues,
                                              parameterLabels,
                                              userSuppliedElements)
      );
      
    }catch(...){
      
      return(arma::datum::nan);
      
    }
    
    return(fitValue);
  }
  
  arma::rowvec gradients(arma::rowvec parameterValues,
                         Rcpp::StringVector parameterLabels) override{
    
    Rcpp::NumericVector gradients(parameterValues.n_elem);
    
    try{
      
      gradients = gradientFunction(parameterValues,
                                   parameterLabels,
                                   userSuppliedElements);
      
    }catch(...){
      
      gradients.fill(NA_REAL);
      return(Rcpp::as<arma::rowvec>(gradients));
      
    }
    
    return(Rcpp::as<arma::rowvec>(gradients)); 
    
  }
  
};

#endif