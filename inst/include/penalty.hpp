#ifndef PENALTY_H
#define PENALTY_H
#include <RcppArmadillo.h>

namespace linr{

template<class T>
class penalty{
public:
  
  virtual double getValue(const Rcpp::NumericVector& parameterValues,
                          const T& tuningParameters);
};
}

#endif