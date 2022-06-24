#ifndef LASSO_H
#define LASSO_H
#include <RcppArmadillo.h>
#include "proximalOperator.hpp"
#include "penalty.hpp"
#include "enet.hpp" // for definition of tuning parameters

namespace linr{

class proximalOperatorLasso: public proximalOperator<tuningParametersEnet>{
public:
  
  arma::rowvec getParameters(const arma::rowvec& parameterValues, 
                                    const arma::rowvec& gradientValues, 
                                    const Rcpp::StringVector& parameterLabels,
                                    const double L,
                                    const tuningParametersEnet& tuningParameters) 
  override {
    
    arma::rowvec u_k = parameterValues - gradientValues/L;
    
    arma::rowvec parameters_kp1(parameterValues.n_elem);
    parameters_kp1.fill(arma::datum::nan);
    
    double lambda_i;
    Rcpp::String parameterLabel;
    int sign;
    for(int p = 0; p < parameterValues.n_elem; p ++)
    {
      
      lambda_i = tuningParameters.alpha *
        tuningParameters.lambda * 
        tuningParameters.weights.at(p);
      
      sign = (u_k.at(p) > 0);
      if(u_k.at(p) < 0) sign = -1;
      parameters_kp1.at(p) = sign*
        std::max(0.0, std::abs(u_k.at(p)) - lambda_i/L);
    }
    return parameters_kp1;
  }
  
};

class penaltyLASSO: public penalty<tuningParametersEnet>{
public:
  
  double getValue(const arma::rowvec& parameterValues, 
                  const Rcpp::StringVector& parameterLabels,
                  const tuningParametersEnet& tuningParameters) 
  override {
    
    double penalty = 0.0;
    double lambda_i;
    
    for(int p = 0; p < parameterValues.n_elem; p ++){
      
      lambda_i = tuningParameters.alpha *
        tuningParameters.lambda * 
        tuningParameters.weights.at(p);
      
      penalty += lambda_i * std::abs(parameterValues.at(p));
    }
    
    return penalty;
  }
  
};

}
#endif
