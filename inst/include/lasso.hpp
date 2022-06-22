#ifndef LASSO_H
#define LASSO_H
#include <RcppArmadillo.h>
#include "proximalOperator.hpp"
#include "penalty.hpp"

namespace linr{

struct tuningParametersLasso{
  double lambda;
  Rcpp::NumericVector weights;
};

class proximalOperatorLasso: public proximalOperator<tuningParametersLasso>{
public:
  
  Rcpp::NumericVector getParameters(const Rcpp::NumericVector& parameterValues, 
                                    const Rcpp::NumericVector& gradientValues, 
                                    const double L,
                                    const tuningParametersLasso& tuningParameters) 
  override {
    
    Rcpp::StringVector parameterLabels = parameterValues.names();
    
    Rcpp::NumericVector u_k = parameterValues - gradientValues/L;
    u_k.names() = parameterLabels;
    Rcpp::NumericVector parameters_kp1(parameterValues.length());
    parameters_kp1.fill(NA_REAL);
    parameters_kp1.names() = parameterLabels;
    
    double lambda_i;
    Rcpp::String parameterLabel;
    int sign;
    for(int p = 0; p < parameterLabels.length(); p ++){
      parameterLabel = parameterLabels.at(p);
      
      lambda_i = tuningParameters.lambda * tuningParameters.weights[parameterLabel];
      
      sign = (u_k[parameterLabel] > 0);
      if(u_k[parameterLabel] < 0) sign = -1;
      parameters_kp1[parameterLabel] = sign*
        std::max(0.0, std::abs(u_k[parameterLabel]) - lambda_i/L);
    }
    return parameters_kp1;
  }
  
};

class penaltyLASSO: public penalty<tuningParametersLasso>{
public:
  
  double getValue(const Rcpp::NumericVector& parameterValues, 
                  const tuningParametersLasso& tuningParameters) 
  override {
    
    Rcpp::StringVector parameterLabels = parameterValues.names();
    
    double penalty = 0.0;
    double lambda_i;
    Rcpp::String parameterLabel;
    
    for(int p = 0; p < parameterLabels.length(); p ++){
      parameterLabel = parameterLabels.at(p);
      
      lambda_i = tuningParameters.lambda * tuningParameters.weights[parameterLabel];
      
      penalty += lambda_i * std::abs(parameterValues[parameterLabel]);
    }
    
    return penalty;
  }
  
};

}
#endif
