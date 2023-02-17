#ifndef LSP_H
#define LSP_H
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

class tuningParametersLSP{
public:
  double lambda;
  double theta;
  arma::rowvec weights;
};

inline double lspPenalty(const double par, 
                         const double lambda,
                         const double theta) {
  
  return(
    lambda * 
      std::log(1.0 + std::abs(par) / theta)
  );
  
}


class proximalOperatorLSP: public proximalOperator<tuningParametersLSP>{
  
public:
  
  arma::rowvec getParameters(const arma::rowvec& parameterValues, 
                             const arma::rowvec& gradientValues, 
                             const Rcpp::StringVector& parameterLabels,
                             const double L,
                             const tuningParametersLSP& tuningParameters) 
  override {
    
    arma::rowvec u_k = parameterValues - gradientValues/L;
    
    arma::rowvec parameters_kp1(parameterValues.n_elem);
    parameters_kp1.fill(arma::datum::nan);
    
    double x, abs_u_k;
    Rcpp::String parameterLabel;
    std::vector<double> C(3, 0.0);
    std::vector<double> xVec(3, 0.0);
    double tempValue;
    int sign;
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++)
    {
      if(tuningParameters.weights.at(p) == 0.0) 
      {
        // unregularized parameter
        parameters_kp1.at(p) = u_k.at(p);
        continue;
      }
      
      abs_u_k = std::abs(u_k.at(p));

      tempValue = std::pow(L,2) * 
        std::pow(abs_u_k - tuningParameters.theta, 2) -
        4.0*L*(tuningParameters.lambda - L*abs_u_k*tuningParameters.theta);
      
      if(tempValue >= 0){
        C.at(1) = std::max(
          (L*(abs_u_k - tuningParameters.theta) + std::sqrt(tempValue)) / (2*L), 
          0.0);
        C.at(2) = std::max(
          (L*(abs_u_k - tuningParameters.theta) - std::sqrt(tempValue)) / (2*L),
          0.0);
        
        for(int c = 0; c < 3; c++){
          
          xVec.at(c) = .5*std::pow(C.at(c) - abs_u_k, 2) +
            (1.0/L)*
            tuningParameters.lambda * 
            std::log(1.0 + C.at(c) / tuningParameters.theta);
          
        }
        
        x = C.at(std::distance(std::begin(xVec), 
                               std::min_element(xVec.begin(), xVec.end())));
        
      }else{
        
        x = 0.0;
        
      }
      
      sign = (u_k.at(p) > 0);
      if(u_k.at(p) < 0) sign = -1;
      
      parameters_kp1.at(p) = sign*x;
      
      
    }
    return parameters_kp1;
  }
  
};

class penaltyLSP: public penalty<tuningParametersLSP>{
public:
  
  double getValue(const arma::rowvec& parameterValues, 
                  const Rcpp::StringVector& parameterLabels,
                  const tuningParametersLSP& tuningParameters) 
  override {
    
    double penalty = 0.0;
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++){
      
      // check if parameter is regularized
      if(tuningParameters.weights.at(p) == 0.0 ) continue;
      
      penalty += lspPenalty(parameterValues.at(p),
                            tuningParameters.lambda,
                            tuningParameters.theta);
    }
    
    return penalty;
  }
  
};

}
#endif
