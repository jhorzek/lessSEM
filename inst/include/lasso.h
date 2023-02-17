#ifndef LASSO_H
#define LASSO_H
#include <RcppArmadillo.h>
#include "proximalOperator.h"
#include "penalty.h"
#include "enet.h" // for definition of tuning parameters

namespace lessSEM{

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
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++)
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
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++){
      
      lambda_i = tuningParameters.alpha *
        tuningParameters.lambda * 
        tuningParameters.weights.at(p);
      
      penalty += lambda_i * std::abs(parameterValues.at(p));
    }
    
    return penalty;
  }
  
  arma::rowvec getSubgradients(const arma::rowvec& parameterValues, 
                               const arma::rowvec& gradients,
                               const tuningParametersEnet& tuningParameters){
    
    arma::rowvec subgradients = gradients;
    double lower, upper;
    int sign;
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p++){
      
      // if not regularized: nothing to do here
      if(tuningParameters.weights.at(p) == 0) continue;
      
      // check if parameter is at non-differentiable place:
      if(parameterValues.at(p) == 0){
        lower = - tuningParameters.weights.at(p) *
          tuningParameters.alpha *
          tuningParameters.lambda;
        // note: we don't add the ridge part here, because this part is already incorporated
        // in the differentiable part in gradients
        upper = -lower;
        
        if(lower < gradients.at(p)){
          subgradients.at(p) = gradients.at(p) + upper;
          continue;
        } else if(gradients.at(p) > upper){
          subgradients.at(p)  = gradients.at(p) + lower;
          continue;
        }else{
          Rcpp::stop("Error in subgradient computation");
        }
        
      }else{
        // parameter is regularized, but not zeroed
        sign = (parameterValues.at(p) > 0);
        if(parameterValues.at(p) < 0) sign = -1;
        
        subgradients.at(p)  = gradients.at(p) +
          sign *
          tuningParameters.weights.at(p) *
          tuningParameters.alpha *
          tuningParameters.lambda;
      }
      
    }// end for parameter
    
    return(subgradients);
    
  }
  
};

}
#endif
