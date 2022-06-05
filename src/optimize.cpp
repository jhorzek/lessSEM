#include <RcppEnsmallen.h>
#include "SEM.hpp"
#include "ensmallen.hpp"
// [[Rcpp::depends(RcppEnsmallen)]]

arma::mat SEMCpp::optimize(arma::mat parameterValues, 
                           const Rcpp::StringVector parameterLabels,
                           const std::string optimizer,
                           double fchange,
                           const bool verbose){
  
  SEMEnsmallen SEMS(*this, parameterLabels, rawData.n_rows, verbose);
  
  if(optimizer.compare("lbfgs") == 0){
    if(verbose){
      Rcpp::Rcout << "Using lbfgs" << std::endl;
    }
    ens::L_BFGS optim;
    optim.Factr() = fchange;
    optim.ArmijoConstant() = .001;
    optim.Optimize(SEMS, parameterValues);
    
  }else if(optimizer.compare("GradientDescent") == 0){
    if(verbose){
      Rcpp::Rcout << "Using GradientDescent" << std::endl;
    }
    
    ens::GradientDescent optim;
    optim.Optimize(SEMS, parameterValues);
    
  }else if(optimizer.compare("DE") == 0){
    if(verbose){
      Rcpp::Rcout << "Using DE" << std::endl;
    }
    ens::DE optim;
    optim.Optimize(SEMS, parameterValues);
    
  }else{
    
    Rcpp::stop("Unknown optimizer selected.");
    
  }
  
  SEMS.clock.stop("timings");
  
  return(parameterValues);
}
