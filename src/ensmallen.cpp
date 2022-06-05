#include <RcppEnsmallen.h>
#include "SEM.hpp"
#include "ensmallen.hpp"
// [[Rcpp::depends(RcppEnsmallen)]]

SEMEnsmallen::SEMEnsmallen(SEMCpp& SEM, const Rcpp::StringVector& parameterLabels): 
  SEM(SEM), parameterLabels(parameterLabels)
{ }

double SEMEnsmallen::Evaluate(const arma::mat& parameterValues){
  
  try{
    SEM.setParameters(parameterLabels, parameterValues, raw);
    SEM.implied();
    SEM.fit();
  }catch(...){
    return(9999999999.99);
  }
  
  if(!SEM.impliedCovariance.is_sympd()){
    return(9999999999.99);
  }
  
  return(SEM.m2LL);
  
}

void SEMEnsmallen::Gradient(const arma::mat& parameterValues, arma::mat& gradients)
{
  try{
    SEM.setParameters(parameterLabels, parameterValues, raw);
    SEM.implied();
    SEM.fit();
    gradients = arma::trans(SEM.getGradients(raw));
  }catch(...){
    gradients.fill(9999999.9999);
    return;
  }
  
  if(!SEM.impliedCovariance.is_sympd()){
    gradients.fill(9999999.9999);
    return;
  }
  
}
