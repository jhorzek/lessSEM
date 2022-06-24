#ifndef SMOOTHPENALTY_H
#define SMOOTHPENALTY_H
#include <RcppArmadillo.h>

namespace linr{

template<class T>
class smoothPenalty{
public:
  
  virtual double getValue(const Rcpp::NumericVector& parameterValues,
                          const T& tuningParameters);
  virtual Rcpp::NumericVector getGradients(const Rcpp::NumericVector& parameterValues,
                                           const T& tuningParameters);
};
}

#endif