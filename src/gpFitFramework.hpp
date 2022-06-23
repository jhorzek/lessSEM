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
  
  double fit(Rcpp::NumericVector parameterValues) override{
    double fitValue;
    try{
      // compute fit
      fitValue = Rcpp::as<double>(fitFunction(parameterValues,
                                              userSuppliedElements)
      );
      
    }catch(...){
      
      return(arma::datum::nan);
      
    }
    
    return(fitValue);
  }
  
  Rcpp::NumericVector gradients(Rcpp::NumericVector parameterValues) override{
    
    Rcpp::NumericVector gradients(parameterValues.length());
    
    try{
      
      gradients = gradientFunction(parameterValues, 
                                   userSuppliedElements);
      
    }catch(...){
      
      gradients.fill(NA_REAL);
      return(gradients);
      
    }
    
    return(gradients); 
    
  }
  
};

#endif