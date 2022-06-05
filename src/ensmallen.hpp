#include <RcppEnsmallen.h>
#include "SEM.hpp"
#include <RcppClock.h>
// [[Rcpp::depends(RcppEnsmallen)]]
//[[Rcpp::depends(RcppClock)]]


class SEMEnsmallen{
public:
  SEMEnsmallen(SEMCpp& SEM, 
               const Rcpp::StringVector& parameterLabels, 
               const double N,
               const bool verbose
  );
  
  double Evaluate(const arma::mat& parameterValues);
  
  void Gradient(const arma::mat& parameterValues, arma::mat& gradients);
  
  Rcpp::Clock clock;
  
private: 
  SEMCpp& SEM;
  const Rcpp::StringVector& parameterLabels;
  const bool raw = true;
  int niter = 0;
  int nInvalid = 0; // count number of invalid iterations
  const bool verbose;
  const double N;
  
};

