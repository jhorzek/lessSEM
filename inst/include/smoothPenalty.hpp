#ifndef SMOOTHPENALTY_H
#define SMOOTHPENALTY_H
#include <RcppArmadillo.h>

namespace linr{

template<class T>
class smoothPenalty{
public:
  
  virtual double getValue(const arma::rowvec& parameterValues,
                          const Rcpp::StringVector& parameterLabels,
                          const T& tuningParameters);
  virtual arma::rowvec getGradients(const arma::rowvec& parameterValues,
                                    const Rcpp::StringVector& parameterLabels,
                                           const T& tuningParameters);
};
}

#endif