#ifndef MODEL_H
#define MODEL_H

#include <RcppArmadillo.h>

namespace ista{

class model{
public:
  virtual double fit(Rcpp::NumericVector parameterValues) = 0;
  virtual Rcpp::NumericVector gradients(Rcpp::NumericVector parameterValues) = 0;
};

}
#endif