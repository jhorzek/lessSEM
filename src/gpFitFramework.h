#ifndef GPFITFRAMEWORK_h
#define GPFITFRAMEWORK_h
#include "SEM.h"
#include "lessSEM.h"

class generalPurposeFitFramework: public lessSEM::model{
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

// creating types for user supplied C++ functions. See
// https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
                Rcpp::List& //additional elements
);
typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;

typedef arma::rowvec (*gradientFunPtr)(const Rcpp::NumericVector&, //parameters
                      Rcpp::List& //additional elements
);
typedef Rcpp::XPtr<gradientFunPtr> gradientFunPtr_t;

class generalPurposeFitFrameworkCpp: public lessSEM::model{
public:
  Rcpp::NumericVector parameters;
  // functions
  fitFunPtr fitFunction;
  gradientFunPtr gradientFunction;
  Rcpp::List userSuppliedElements;
  
  generalPurposeFitFrameworkCpp(
    Rcpp::NumericVector parameters_,
    SEXP fitFunctionSEXP,
    SEXP gradientFunctionSEXP,
    Rcpp::List userSuppliedElements_
  ){
    
    parameters = parameters_;
    
    Rcpp::XPtr<fitFunPtr> xpFitFun(fitFunctionSEXP);
    fitFunction = *xpFitFun;
    
    Rcpp::XPtr<gradientFunPtr> xpGradientFun(gradientFunctionSEXP);
    gradientFunction = *xpGradientFun;
    
    userSuppliedElements = userSuppliedElements_;
    
  }
  
  double fit(arma::rowvec parameterValues,
             Rcpp::StringVector parameterLabels) override{
               for(int i = 0; i < parameterValues.n_elem; i++){
                 parameters.at(i) = parameterValues.at(i);
               }
               
               double fitValue;
               try{
                 // compute fit
                 fitValue = fitFunction(parameters,
                                        userSuppliedElements);
                 
               }catch(...){
                 
                 return(arma::datum::nan);
                 
               }
               
               return(fitValue);
             }
  
  arma::rowvec gradients(arma::rowvec parameterValues,
                         Rcpp::StringVector parameterLabels) override{
                           
                           for(int i = 0; i < parameterValues.n_elem; i++){
                             parameters.at(i) = parameterValues.at(i);
                           }
                           
                           Rcpp::NumericVector gradients(parameterValues.n_elem);
                           
                           try{
                             
                             gradients = gradientFunction(parameters,
                                                          userSuppliedElements);
                             
                           }catch(...){
                             
                             gradients.fill(NA_REAL);
                             return(Rcpp::as<arma::rowvec>(gradients));
                             
                           }
                           
                           return(Rcpp::as<arma::rowvec>(gradients)); 
                           
                         }
  
};

#endif