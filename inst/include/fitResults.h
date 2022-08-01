#ifndef FITRESULTS_H
#define FITRESULTS_H

namespace lessSEM{

struct fitResults{
  double fit; // the final fit value (regularized fit)
  arma::rowvec fits; // a vector with all fits at the outer iteration
  bool convergence; // was the outer breaking condition met?
  arma::rowvec parameterValues; // final parameter values
  arma::mat Hessian;
};

}

#endif