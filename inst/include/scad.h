#ifndef SCAD2_H
#define SCAD2_H
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

class tuningParametersScad{
public:
  double lambda;
  double theta;
  arma::rowvec weights;
};

inline double scadPenalty(const double par, 
                           const double lambda,
                           const double theta){
  
  double absPar = std::abs(par);
  
  if(absPar <= lambda){
    // reduces to lasso penalty
    return(lambda * absPar);
    
  }else if((lambda < absPar) && (absPar <= lambda * theta)){
    // reduces to a smooth penalty
    return(
      (-std::pow(par, 2) + 
        2.0 * theta * lambda * absPar - std::pow(lambda,2))/
          (2.0 * (theta - 1.0))
    );
    
  }else if( absPar >  (lambda * theta)){
    // reduces to a constant penalty
    return( (( theta + 1.0 ) * std::pow(lambda,2)) / 2.0);
    
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
    
    double abs_u_k, v;
    Rcpp::String parameterLabel;
    std::vector<double> x(4, 0.0); // to save the minima of the 
    // three different regions of the penalty function
    std::vector<double> h(4, 0.0); // to save the function values of the 
    // four possible minima saved in x
    int sign; // sign of u_k
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
      
      abs_u_k = std::abs(u_k.at(p));
      
      // assume that the solution is found in 
      // |x| <= lambda. In this region, the 
      // scad penalty is identical to the lasso, so
      // we can use the same minimizer as for the lasso
      // with the additional bound that 
      
      // identical to Gong et al. (2013)
      x.at(0) = sign * std::min(
        tuningParameters.lambda, 
        std::max(
          0.0,
          abs_u_k - tuningParameters.lambda/L
        )
      );
      
      // assume that lambda <= |u| <= theta*lambda
      // The following differs from Gong et al. (2013)
      
      v = 1.0 - 1.0/(L*(tuningParameters.theta - 1.0)); // used repeatedly;
      // only computed for convenience
      
      x.at(1) = std::min(
        thetaXlambda, 
        std::max(
          tuningParameters.lambda,
          (u_k.at(p) / v) - (thetaXlambda)/(L*(tuningParameters.theta - 1.0)*v)
        )
      );
      
      x.at(2) = std::max(
        -thetaXlambda, 
        std::min(
          -tuningParameters.lambda,
          (u_k.at(p) / v) + (thetaXlambda)/(L*(tuningParameters.theta - 1.0)*v)
        )
      );
      
      // assume that |u| >= lambda*theta
      // identical to Gong et al. (2013)
      
      x.at(3) = sign * std::max(
        thetaXlambda,
        abs_u_k
      );
      
      for(int i = 0; i < 4; i++){
        h.at(i) = .5 * std::pow(x.at(i) - u_k.at(p), 2) + // distance between parameters
          (1.0/L) * scadPenalty(x.at(i), tuningParameters.lambda, tuningParameters.theta);
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
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++){
      // unregularized values:
      if(tuningParameters.weights.at(p) == 0.0) continue;
      
      // scad penalty value:
      penalty += scadPenalty(parameterValues.at(p), 
                              tuningParameters.lambda,
                              tuningParameters.theta);
      
    }
    
    return penalty;
  }
  
};

}
#endif
