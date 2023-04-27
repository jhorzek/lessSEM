#ifndef SEMFF_h
#define SEMFF_h
#include "SEM.h"
#include "lessSEM.h"

template<typename sem>
class SEMFitFramework: public lessSEM::model{
public:
  
  sem& SEM;
  const double scaleFit; // option to up - or downscale the fit function values
  
  SEMFitFramework(sem& SEM_,
                  double scaleFit_ = 1.0): SEM(SEM_), scaleFit(scaleFit_){}
  
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
               if(!SEM.impliedIsPD()){
                 return(arma::datum::nan);
               }
               return(scaleFit * SEM.m2LL);
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
                           if(!SEM.impliedIsPD()){
                             gradients.fill(arma::datum::nan);
                             return(gradients);
                           }
                           
                           return(scaleFit * gradients); 
                           
                         }
  
};

#endif