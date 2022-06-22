#ifndef PROXIMALOPERATOR_H
#define PROXIMALOPERATOR_H
#include <RcppArmadillo.h>

namespace linr{
template<class T>
class proximalOperator{
public:
  virtual Rcpp::NumericVector getParameters(const Rcpp::NumericVector& parameterValues, 
                                            const Rcpp::NumericVector& gradientValues,
                                            const double L,
                                            const T& tuningParameters) = 0;
};
}
#endif