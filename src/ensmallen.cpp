#include <RcppEnsmallen.h>
#include "SEM.hpp"
#include "ensmallen.hpp"
// [[Rcpp::depends(RcppEnsmallen)]]

SEMEnsmallen::SEMEnsmallen(SEMCpp& SEM, 
                           const Rcpp::StringVector& parameterLabels, 
                           const double N,
                           const bool verbose): 
  SEM(SEM), parameterLabels(parameterLabels), N(N), verbose(verbose)
{ }

double SEMEnsmallen::Evaluate(const arma::mat& parameterValues){
  niter++;
  try{
    clock.tick("parameters");
    SEM.setParameters(parameterLabels, parameterValues, raw);
    clock.tock("parameters");
    clock.tick("implied");
    SEM.implied();
    clock.tock("implied");
    clock.tick("fit");
    SEM.fit();
    clock.tock("fit");
  }catch(...){
    nInvalid++;
    return(9999999999.99);
  }
  
  if(!SEM.impliedCovariance.is_sympd()){
    nInvalid++;
    return(9999999999.99);
  }
  if(verbose && niter%10 == 0){
    Rcpp::Rcout << "Iteration " << niter << ": " << SEM.m2LL << " (" << nInvalid << " invalid iterations)" << std::endl;
  }
  return(SEM.m2LL/N);
  
}

void SEMEnsmallen::Gradient(const arma::mat& parameterValues, arma::mat& gradients)
{
  try{
    clock.tick("parameters_grd");
    SEM.setParameters(parameterLabels, parameterValues, raw);
    clock.tock("parameters_grd");
    clock.tick("implied_grd");
    SEM.implied();
    clock.tock("implied_grd");
    clock.tick("fit_grd");
    SEM.fit();
    clock.tock("fit_grd");
    clock.tick("gradients_grd");
    gradients = arma::trans(SEM.getGradients(raw)/N);
    clock.tock("gradients_grd");
  }catch(...){
    nInvalid++;
    gradients.fill(9999999.9999);
    return;
  }
  
  if(!SEM.impliedCovariance.is_sympd()){
    gradients.fill(9999999.9999);
    return;
  }
  
}
