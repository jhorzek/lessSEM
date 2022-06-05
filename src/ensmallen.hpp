#include <RcppEnsmallen.h>
#include "SEM.hpp"
// [[Rcpp::depends(RcppEnsmallen)]]

class SEMEnsmallen{
  public:
    SEMEnsmallen(SEMCpp& SEM, const Rcpp::StringVector& parameterLabels);
  
  double Evaluate(const arma::mat& parameterValues);
  
  void Gradient(const arma::mat& parameterValues, arma::mat& gradients);
  
  private: 
    SEMCpp& SEM;
  const Rcpp::StringVector& parameterLabels;
  const bool raw = true;
  
};

