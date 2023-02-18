#ifndef CAPPEDL1_H
#define CAPPEDL1_H
#include <RcppArmadillo.h>
#include "proximalOperator.h"
#include "penalty.h"
#include "enet.h" // for definition of tuning parameters

// The proximal operator for this penalty function has been developed by
// Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
// A general iterative shrinkage and thresholding algorithm for non-convex 
// regularized optimization problems. Proceedings of the 30th International 
// Conference on Machine Learning, 28(2)(2), 37–45.
// The implementation directly follows that of Gong et al. (2013)

namespace lessSEM{
  
  class tuningParametersCappedL1: public lessSEM::tuningParametersEnet{
  public:
    double theta; // threshold parameter; any parameter
    // above this threshold will only receive the constant penalty
    // lambda_i*theta, all below will get lambda_i*parameterValue_i
  };
  
  class proximalOperatorCappedL1: public proximalOperator<tuningParametersCappedL1>{
    
    public:
      
      arma::rowvec getParameters(const arma::rowvec& parameterValues, 
                                 const arma::rowvec& gradientValues, 
                                 const Rcpp::StringVector& parameterLabels,
                                 const double L,
                                 const tuningParametersCappedL1& tuningParameters) 
    override {
      
      // step in descending direction with step size (1/L):
      arma::rowvec u_k = parameterValues - gradientValues/L;
      
      arma::rowvec parameters_kp1(parameterValues.n_elem);
      parameters_kp1.fill(arma::datum::nan);
      
      double lambda_i, x_1, x_2, h_1, h_2, abs_u_k;
      Rcpp::String parameterLabel;
      int sign;
      for(unsigned int p = 0; p < parameterValues.n_elem; p ++)
      {
        
        lambda_i = tuningParameters.alpha *
          tuningParameters.lambda * 
          tuningParameters.weights.at(p);
        
        sign = (u_k.at(p) > 0);
        if(u_k.at(p) < 0) sign = -1;
        
        abs_u_k = std::abs(u_k.at(p));
        
        x_1 = sign*std::max(abs_u_k, tuningParameters.theta);
        x_2 = sign*std::min(tuningParameters.theta,
                            std::max(abs_u_k - lambda_i/L, 0.0));
        // h_1 and h_2 will always be positive. The minimum is therefore 
        // 0 which is also the value we get if either x_1 or x_2 are 
        // equivalent to the proposed parameter u_k in descend-direction.
        // This is the case if the absolute value of the
        // proposed parameter is above the threshold theta -> x_1 = u_k.
        // => IF |u_k| > THETA, WE ALWAYS SELECT u_k
        // If the proposed parameter |u_k| is below the threshold theta
        // x_2 comes into play. x_2 is at minimum equal to theta (upper bound)
        // and otherwise equal to std::max(abs_u_k - lambda_i/L, 0.0)
        // which is the proximal operator of the lasso penalty
        // => IF |u_k| > THETA, WE ALWAYS TAKE THE NORMAL LASSO UPDATE
        h_1 = .5*std::pow(x_1 - u_k.at(p), 2) + 
          (lambda_i/L)*std::min(std::abs(x_1), tuningParameters.theta);
        h_2 = .5*std::pow(x_2 - u_k.at(p), 2) + 
          (lambda_i/L)*std::min(std::abs(x_2), tuningParameters.theta);
        
        if(h_1 <= h_2){
          parameters_kp1.at(p) = x_1;
        }else{
          parameters_kp1.at(p) = x_2;
        }
      }
      return parameters_kp1;
    }
    
  };
  
  class penaltyCappedL1: public penalty<tuningParametersCappedL1>{
    public:
      
      double getValue(const arma::rowvec& parameterValues, 
                      const Rcpp::StringVector& parameterLabels,
                      const tuningParametersCappedL1& tuningParameters) 
    override {
      
      double penalty = 0.0;
      double lambda_i;
      
      for(unsigned int p = 0; p < parameterValues.n_elem; p ++){
        
        lambda_i = tuningParameters.alpha *
          tuningParameters.lambda * 
          tuningParameters.weights.at(p);
        
        penalty += lambda_i * std::min(std::abs(parameterValues.at(p)), 
                                       tuningParameters.theta);
      }
      
      return penalty;
    }
  
};

}
#endif
