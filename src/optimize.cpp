#include <RcppEnsmallen.h>
#include "SEM.hpp"
#include "ensmallen.hpp"
// [[Rcpp::depends(RcppEnsmallen)]]

arma::mat SEMCpp::optimize(arma::mat parameterValues, 
                           const Rcpp::StringVector parameterLabels,
                           const std::string optimizer){
  
  SEMEnsmallen SEMS(*this, parameterLabels);
  
  if(optimizer.compare("lbfgs") == 0){
    
    ens::L_BFGS optim;
    optim.Optimize(SEMS, parameterValues);
    
  }else{
    
    Rcpp::stop("Unknown optimizer selected.");
    
  }
  
  return(parameterValues);
}
