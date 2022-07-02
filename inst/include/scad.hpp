#ifndef SCAD_H
#define SCAD_H
#include <RcppArmadillo.h>
#include "proximalOperator.hpp"
#include "penalty.hpp"

// The proximal operator for this penalty function has been developed by
// Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
// A general iterative shrinkage and thresholding algorithm for non-convex 
// regularized optimization problems. Proceedings of the 30th International 
// Conference on Machine Learning, 28(2)(2), 37–45.
// The implementation directly follows that of Gong et al. (2013)

namespace lessSEM{

class tuningParametersScad{
public:
  double lambda;
  double theta;
  arma::rowvec weights;
};

inline double scadPenalty(const double par, 
                          const double lambda_p,
                          const double theta){
  
  double absPar = std::abs(par);
  
  if(absPar <= lambda_p){
    
    return(lambda_p * absPar);
    
  }else if((lambda_p < absPar) & (absPar <= lambda_p * theta)){
    
    return(
      (-std::pow(par, 2) + 
        2.0 * theta * lambda_p * absPar - std::pow(lambda_p,2))/
          (2.0 * (theta - 1.0))
    );
    
  }else if( absPar >  (lambda_p * theta)){
    
    return( ( theta + 1 ) * std::pow(lambda_p,2) / 2.0);
    
  }else{
    
    Rcpp::stop("Error while evaluating scad");
    
  }
  
  // the following should never be called:
  return 0.0;
  
}


class proximalOperatorScad: public proximalOperator<tuningParametersScad>{
  
public:
  
  arma::rowvec getParameters(const arma::rowvec& parameterValues, 
                             const arma::rowvec& gradientValues, 
                             const Rcpp::StringVector& parameterLabels,
                             const double L,
                             const tuningParametersScad& tuningParameters) 
  override {
    
    arma::rowvec u_k = parameterValues - gradientValues/L;
    
    arma::rowvec parameters_kp1(parameterValues.n_elem);
    parameters_kp1.fill(arma::datum::nan);
    
    double lambda_p, theta_p, abs_u_k;
    Rcpp::String parameterLabel;
    std::vector<double> x(3, 0.0);
    std::vector<double> h(3, 0.0);
    int sign;
    
    for(int p = 0; p < parameterValues.n_elem; p ++)
    {
      
      lambda_p = tuningParameters.weights.at(p) * tuningParameters.lambda;
      theta_p = tuningParameters.weights.at(p) * tuningParameters.theta;
      
      sign = (u_k.at(p) > 0);
      if(u_k.at(p) < 0) sign = -1;
      
      abs_u_k = std::abs(u_k.at(p));
      
      x.at(0) = sign*std::min(lambda_p, std::max(
        0.0,
        abs_u_k - lambda_p/L
      ));
      
      x.at(1) = sign*std::min(theta_p * lambda_p, 
           std::max(
             lambda_p,
             (L * abs_u_k*(theta_p - 1) - 
               theta_p*lambda_p)/
                 ( L*(theta_p - 2.0) )
           )
      );
      
      x.at(2) = sign*std::max(theta_p * lambda_p, 
           abs_u_k);
      
      double penalty;
      for(int c = 0; c < 3; c++){
        
        // scad penalty value:
        penalty = scadPenalty(
          x.at(c), 
          lambda_p,
          theta_p
        );
        
        h.at(c) = .5*std::pow(x.at(c) - u_k.at(p), 2) // distance between parameters
          + (1.0/L) * penalty; // add penalty value of scad:
        
      }
      parameters_kp1.at(p) = x.at(std::distance(std::begin(h), 
                                  std::min_element(h.begin(), h.end())));
      
    }
    
    return parameters_kp1;
  }
  
};

class penaltyScad: public penalty<tuningParametersScad>{
public:
  
  double getValue(const arma::rowvec& parameterValues, 
                  const Rcpp::StringVector& parameterLabels,
                  const tuningParametersScad& tuningParameters) 
  override {
    
    double penalty = 0.0;
    double absPar, lambda_p, theta_p;
    for(int p = 0; p < parameterValues.n_elem; p ++){
      
      lambda_p = tuningParameters.weights.at(p) * tuningParameters.lambda;
      theta_p = tuningParameters.weights.at(p) * tuningParameters.theta;
      
      // scad penalty value:
      penalty += scadPenalty(parameterValues.at(p), 
                             lambda_p,
                             tuningParameters.theta);
      
    }
    
    return penalty;
  }
  
};

}
#endif
