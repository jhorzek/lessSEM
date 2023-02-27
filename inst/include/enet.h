#ifndef ENET_H
#define ENET_H
#include <RcppArmadillo.h>

// The elastic net is a combination of ridge and lasso
// In this file, only the tuning parameters required by lasso
// and ridge are defined.
namespace lessSEM{
  
  class tuningParametersEnet{
  public:
    double lambda;
    double alpha;
    arma::rowvec weights;
  };
  
  // for glmnet, we will allow for different alphas and lambdas to combine penalties
  class tuningParametersEnetGlmnet{
  public:
    arma::rowvec lambda;
    arma::rowvec alpha;
    arma::rowvec weights;
  };
}

#endif