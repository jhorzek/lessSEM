#ifndef MIXEDPENALTY_H
#define MIXEDPENALTY_H
#include <RcppArmadillo.h>
#include "proximalOperator.h"
#include "penalty.h"
#include "cappedL1.h"
#include "lasso.h"
#include "lsp.h"
#include "mcp.h"
#include "scad.h"


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

enum penaltyType{
  none,
  cappedL1,
  lasso,
  lsp,
  mcp,
  scad
};
const std::vector<std::string> penaltyType_txt = {
  "none",
  "cappedL1",
  "lasso",
  "lsp",
  "mcp",
  "scad"
};

class tuningParametersMixedPenalty{
public:
  arma::rowvec lambda;
  arma::rowvec theta;
  arma::rowvec alpha;
  arma::rowvec weights;
  std::vector<penaltyType> pt;
};

class proximalOperatorMixedPenalty: public proximalOperator<tuningParametersMixedPenalty>{
  
public:
  
  arma::rowvec getParameters(const arma::rowvec& parameterValues, 
                             const arma::rowvec& gradientValues, 
                             const Rcpp::StringVector& parameterLabels,
                             const double L,
                             const tuningParametersMixedPenalty& tuningParameters) 
  override {
    
    arma::rowvec u_k = parameterValues - gradientValues/L;
    
    arma::rowvec parameters_kp1(parameterValues.n_elem);
    parameters_kp1.fill(arma::datum::nan);
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++)
    {
      
      if(tuningParameters.weights.at(p) == 0.0) 
      {
        // unregularized parameter
        parameters_kp1.at(p) = u_k.at(p);
        continue;
      }
      
      if(tuningParameters.pt.at(p) == cappedL1){
        
        double lambda_i, x_1, x_2, h_1, h_2, abs_u_k;
        Rcpp::String parameterLabel;
        int sign;
        lambda_i = tuningParameters.alpha.at(p) *
          tuningParameters.lambda.at(p) * 
          tuningParameters.weights.at(p);
        
        sign = (u_k.at(p) > 0);
        if(u_k.at(p) < 0) sign = -1;
        
        abs_u_k = std::abs(u_k.at(p));
        
        x_1 = sign*std::max(abs_u_k, tuningParameters.theta.at(p));
        x_2 = sign*std::min(tuningParameters.theta.at(p),
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
          (lambda_i/L)*std::min(std::abs(x_1), tuningParameters.theta.at(p));
        h_2 = .5*std::pow(x_2 - u_k.at(p), 2) + 
          (lambda_i/L)*std::min(std::abs(x_2), tuningParameters.theta.at(p));
        
        if(h_1 <= h_2){
          parameters_kp1.at(p) = x_1;
        }else{
          parameters_kp1.at(p) = x_2;
        }
        
        continue;
      }
      
      
      if(tuningParameters.pt.at(p) == lasso){
        
        double lambda_i;
        Rcpp::String parameterLabel;
        int sign;
        
        lambda_i = tuningParameters.lambda.at(p) * 
          tuningParameters.weights.at(p);
        
        sign = (u_k.at(p) > 0);
        if(u_k.at(p) < 0) sign = -1;
        parameters_kp1.at(p) = sign*
          std::max(0.0, std::abs(u_k.at(p)) - lambda_i/L);
        
        continue;
      }
      
      if(tuningParameters.pt.at(p) == lsp){
        
        double x, abs_u_k;
        Rcpp::String parameterLabel;
        std::vector<double> C(3, 0.0);
        std::vector<double> xVec(3, 0.0);
        double tempValue;
        int sign;
        
        abs_u_k = std::abs(u_k.at(p));
        
        tempValue = std::pow(L,2) * 
          std::pow(abs_u_k - tuningParameters.theta.at(p), 2) -
          4.0*L*(tuningParameters.lambda.at(p) - L*abs_u_k*tuningParameters.theta.at(p));
        
        if(tempValue >= 0){
          C.at(1) = std::max(
            (L*(abs_u_k - tuningParameters.theta.at(p)) + std::sqrt(tempValue)) / (2*L), 
            0.0);
          C.at(2) = std::max(
            (L*(abs_u_k - tuningParameters.theta.at(p)) - std::sqrt(tempValue)) / (2*L),
            0.0);
          
          for(int c = 0; c < 3; c++){
            
            xVec.at(c) = .5*std::pow(C.at(c) - abs_u_k, 2) +
              (1.0/L)*
              tuningParameters.lambda.at(p) * 
              std::log(1.0 + C.at(c) / tuningParameters.theta.at(p));
            
          }
          
          x = C.at(std::distance(std::begin(xVec), 
                                 std::min_element(xVec.begin(), xVec.end())));
          
        }else{
          
          x = 0.0;
          
        }
        
        sign = (u_k.at(p) > 0);
        if(u_k.at(p) < 0) sign = -1;
        
        parameters_kp1.at(p) = sign*x;
        
        continue;
      }
      
      if(tuningParameters.pt.at(p) == mcp){
        
        double abs_u_k, v;
        std::vector<double> x(4, 0.0);
        std::vector<double> h(4, 0.0);
        int sign;
        double thetaXlambda = tuningParameters.theta.at(p)*tuningParameters.lambda.at(p);
        
        sign = (u_k.at(p) > 0);
        if(u_k.at(p) < 0) sign = -1;
        
        v = 1.0 - 1.0/(L*tuningParameters.theta.at(p)); // used repeatedly;
        // only computed for convenience
        
        abs_u_k = std::abs(u_k.at(p));
        
        // Assume that x = 0
        x.at(0) = 0.0;
        
        // Assume that x > 0 and x <= theta*lambda
        x.at(1) = std::min(
          thetaXlambda,
          u_k.at(p)/v - 1.0/(L*v)*tuningParameters.lambda.at(p)
        );
        
        // Assume that x < 0 and x => - theta*lambda
        x.at(2) = std::max(
          -thetaXlambda,
          u_k.at(p)/v + 1.0/(L*v)*tuningParameters.lambda.at(p)
        );
        
        // Assume that |x| >  theta*lambda
        x.at(3) = sign*std::max(
          thetaXlambda,
          abs_u_k
        );
        
        for(int i = 0; i < 4; i++){
          h.at(i) = .5 * std::pow(x.at(i) - u_k.at(p), 2) + // distance between parameters
            (1.0/L) * mcpPenalty(x.at(i), 
             tuningParameters.lambda.at(p), 
             tuningParameters.theta.at(p));
        }
        
        parameters_kp1.at(p) = x.at(std::distance(std::begin(h), 
                                    std::min_element(h.begin(), h.end())));
        
        continue;
      }
      
      if(tuningParameters.pt.at(p) == scad){
        
        double abs_u_k, v;
        Rcpp::String parameterLabel;
        std::vector<double> x(4, 0.0); // to save the minima of the 
        // three different regions of the penalty function
        std::vector<double> h(4, 0.0); // to save the function values of the 
        // four possible minima saved in x
        int sign; // sign of u_k
        double thetaXlambda = tuningParameters.theta.at(p)*tuningParameters.lambda.at(p);
        
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
          tuningParameters.lambda.at(p), 
          std::max(
            0.0,
            abs_u_k - tuningParameters.lambda.at(p)/L
          )
        );
        
        // assume that lambda <= |u| <= theta*lambda
        // The following differs from Gong et al. (2013)
        
        v = 1.0 - 1.0/(L*(tuningParameters.theta.at(p) - 1.0)); // used repeatedly;
        // only computed for convenience
        
        x.at(1) = std::min(
          thetaXlambda, 
          std::max(
            tuningParameters.lambda.at(p),
            (u_k.at(p) / v) - (thetaXlambda)/(L*(tuningParameters.theta.at(p) - 1.0)*v)
          )
        );
        
        x.at(2) = std::max(
          -thetaXlambda, 
          std::min(
            -tuningParameters.lambda.at(p),
            (u_k.at(p) / v) + (thetaXlambda)/(L*(tuningParameters.theta.at(p) - 1.0)*v)
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
            (1.0/L) * scadPenalty(x.at(i), tuningParameters.lambda.at(p), tuningParameters.theta.at(p));
        }
        
        parameters_kp1.at(p) = x.at(std::distance(std::begin(h), 
                                    std::min_element(h.begin(), h.end())));
        
        continue;
      }
      
      Rcpp::stop("Unknown penalty function.");
      
      
    }
    
    return parameters_kp1;
  }
  
};

class penaltyMixedPenalty: public penalty<tuningParametersMixedPenalty>{
public:
  
  double getValue(const arma::rowvec& parameterValues, 
                  const Rcpp::StringVector& parameterLabels,
                  const tuningParametersMixedPenalty& tuningParameters) 
  override {
    
    double penalty = 0.0;
    double lambda_i;
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++){
      // unregularized values:
      if(tuningParameters.weights.at(p) == 0.0) continue;
      
      switch(tuningParameters.pt.at(p)){
      
      case cappedL1:
        lambda_i = tuningParameters.alpha.at(p) *
          tuningParameters.lambda.at(p) * 
          tuningParameters.weights.at(p);
        
        penalty += lambda_i * std::min(std::abs(parameterValues.at(p)), 
                                       tuningParameters.theta.at(p));
        
        break;
        
      case lasso:
        
        lambda_i = tuningParameters.lambda.at(p) * 
          tuningParameters.weights.at(p);
        
        penalty += lambda_i * std::abs(parameterValues.at(p));
        break;
        
      case lsp:
        
        penalty += lspPenalty(parameterValues.at(p),
                              tuningParameters.lambda.at(p),
                              tuningParameters.theta.at(p));
        
        break;
        
      case mcp: 
        
        penalty += mcpPenalty(parameterValues.at(p), 
                              tuningParameters.lambda.at(p),
                              tuningParameters.theta.at(p));
        
        break;
        
      case scad:  
        
        penalty += scadPenalty(parameterValues.at(p), 
                               tuningParameters.lambda.at(p),
                               tuningParameters.theta.at(p));
        
        break;
        
      default:
        Rcpp::stop("Penalty not found.");
      }
      
    }
    
    return penalty;
  }
  
};

}
#endif
