#ifndef MODEL_H
#define MODEL_H

#include <RcppArmadillo.h>

namespace lessSEM{

class model{
public:
  virtual double fit(arma::rowvec parameterValues,
                     Rcpp::StringVector parameterLabels) = 0;
  virtual arma::rowvec gradients(arma::rowvec parameterValues, 
                                 Rcpp::StringVector parameterLabels) = 0;
};

}
#endif