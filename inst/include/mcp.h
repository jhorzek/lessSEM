#ifndef MCP_H
#define MCP_H
#include <RcppArmadillo.h>
#include "proximalOperator.h"
#include "penalty.h"

// The proximal operator for this penalty function has been developed by
// Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
// A general iterative shrinkage and thresholding algorithm for non-convex 
// regularized optimization problems. Proceedings of the 30th International 
// Conference on Machine Learning, 28(2)(2), 37–45.

// The implementation deviated from that of Gong et al. (2013). We 
// could not get the procedure outlined by Gong et al. (2013) to return the 
// correct results. In the following, we will use our own derivation, but
// it is most likely equivalent to that of Gong et al. (2013) and we 
// just failed to implement their approach correctly. If you find our
// mistake, feel free to improve upon the following implementation.

namespace lessSEM{

class tuningParametersMcp{
public:
  double lambda;
  double theta;
  arma::rowvec weights;
};

inline double mcpPenalty(const double par, 
                         const double lambda,
                         const double theta) {
  
  double absPar = std::abs(par);
  
  if(absPar <= (lambda * theta)){
    return(
      lambda * absPar - std::pow( par, 2 ) / (2.0 * theta)
    );
    
  }else if(absPar > (lambda * theta)){
    
    return(
      theta * std::pow(lambda, 2) / 2.0
    );
    
  }else{
    Rcpp::stop("Error while evaluating mcp");
  }
  // The following should never be called:
  return 0.0;
}

class proximalOperatorMcp: public proximalOperator<tuningParametersMcp>{
  
public:
  
  arma::rowvec getParameters(const arma::rowvec& parameterValues, 
                             const arma::rowvec& gradientValues, 
                             const Rcpp::StringVector& parameterLabels,
                             const double L,
                             const tuningParametersMcp& tuningParameters) 
  override {
    
    arma::rowvec u_k = parameterValues - gradientValues/L;
    
    arma::rowvec parameters_kp1(parameterValues.n_elem);
    parameters_kp1.fill(arma::datum::nan);
    
    double abs_u_k, v;
    std::vector<double> x(4, 0.0);
    std::vector<double> h(4, 0.0);
    int sign;
    double thetaXlambda = tuningParameters.theta*tuningParameters.lambda;
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++)
    {
      
      if(tuningParameters.weights.at(p) == 0.0) 
      {
        // unregularized parameter
        parameters_kp1.at(p) = u_k.at(p);
        continue;
      }
      
      sign = (u_k.at(p) > 0);
      if(u_k.at(p) < 0) sign = -1;
      
      v = 1.0 - 1.0/(L*tuningParameters.theta); // used repeatedly;
      // only computed for convenience
      
      abs_u_k = std::abs(u_k.at(p));
      
      // Assume that x = 0
      x.at(0) = 0.0;
      
      // Assume that x > 0 and x <= theta*lambda
      x.at(1) = std::min(
        thetaXlambda,
        u_k.at(p)/v - 1.0/(L*v)*tuningParameters.lambda
      );
      
      // Assume that x < 0 and x => - theta*lambda
      x.at(2) = std::max(
        -thetaXlambda,
        u_k.at(p)/v + 1.0/(L*v)*tuningParameters.lambda
      );
      
      // Assume that |x| >  theta*lambda
      x.at(3) = sign*std::max(
        thetaXlambda,
        abs_u_k
      );
      
      for(int i = 0; i < 4; i++){
        h.at(i) = .5 * std::pow(x.at(i) - u_k.at(p), 2) + // distance between parameters
          (1.0/L) * mcpPenalty(x.at(i), 
           tuningParameters.lambda, 
           tuningParameters.theta);
      }
      
      parameters_kp1.at(p) = x.at(std::distance(std::begin(h), 
                                  std::min_element(h.begin(), h.end())));
      
    }
    return parameters_kp1;
  }
  
};

class penaltyMcp: public penalty<tuningParametersMcp>{
public:
  
  double getValue(const arma::rowvec& parameterValues, 
                  const Rcpp::StringVector& parameterLabels,
                  const tuningParametersMcp& tuningParameters) 
  override {
    
    double penalty = 0.0;
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++){
      
      // unregularized values:
      if(tuningParameters.weights.at(p) == 0.0) continue;
      
      // mcp penalty value:
      penalty += mcpPenalty(parameterValues.at(p), 
                             tuningParameters.lambda,
                             tuningParameters.theta);
      
    }
    
    return penalty;
  }
  
};

}
#endif
