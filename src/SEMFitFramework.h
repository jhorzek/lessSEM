#ifndef SEMFF_h
#define SEMFF_h
#include "SEM.h"
#include "lessSEM.h"

template<typename sem>
class SEMFitFramework: public lessSEM::model{
public:
  
  sem& SEM;
  
  SEMFitFramework(sem& SEM_): SEM(SEM_){}
  
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
                           if(!SEM.impliedIsPD()){
                             gradients.fill(arma::datum::nan);
                             return(gradients);
                           }
                           
                           return(gradients); 
                           
                         }
  
};

#endif