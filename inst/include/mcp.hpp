#ifndef MCP_H
#define MCP_H
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
  
  class tuningParametersMcp{
    public:
      double lambda;
    double theta;
    arma::rowvec weights;
  };
  
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
      
      double lambda_i, abs_u_k, h_1, h_2;
      Rcpp::String parameterLabel;
      std::vector<double> C(3, 0.0);
      std::vector<double> f(3, 0.0);
      std::vector<double> x(2, 0.0);
      int sign, nval;
      
      for(int p = 0; p < parameterValues.n_elem; p ++)
      {
        
        lambda_i = tuningParameters.weights.at(p) * tuningParameters.lambda;
        
        sign = (u_k.at(p) > 0);
        if(u_k.at(p) < 0) sign = -1;
        
        abs_u_k = std::abs(u_k.at(p));
        
        // compute z first; needed for x_1
        C.at(1) = lambda_i * tuningParameters.theta;
        
        if(tuningParameters.theta - 1 != 0){
        
        C.at(2) = std::min(
          lambda_i * tuningParameters.theta,
          std::max(
            0.0,
            (tuningParameters.theta *(L*abs_u_k - lambda_i)) / 
              (L * (tuningParameters.theta - 1.0))
          )
        );
          
          nval = 3;
          
        }else{
          
          // only the first two values are valid otherwise
          nval = 2;
          
        }
        
        for(int c = 0; c < nval; c++){
          f.at(c) = .5*std::pow( C.at(c) - abs_u_k , 2) + 
            (lambda_i/L) * C.at(c) - std::pow(C.at(c),2)/(2*tuningParameters.theta);
        }
        
        x.at(0) = sign * C.at(std::distance(std::begin(f), 
                           std::min_element(f.begin(), f.end() - (3-nval))));
        
        x.at(1) = sign * std::max( tuningParameters.theta * lambda_i, abs_u_k );
        
        double absPar, penalty;
        for(int c = 0; c < 2; c++){
          
          absPar = std::abs(x.at(c));
          
          if(absPar <= lambda_i * tuningParameters.theta){
            
            penalty = lambda_i * absPar - std::pow( x.at(c), 2 ) / (2*tuningParameters.theta);
            
          }else if(absPar > lambda_i * tuningParameters.theta){
            
            penalty = tuningParameters.theta * std::pow(lambda_i, 2) / 2.0;
            
          }else{
            Rcpp::stop("Error while evaluating mcp");
          }
          
          f.at(c) = .5*std::pow( x.at(c) - abs_u_k , 2) + 
            (lambda_i/L) * penalty;
        }
        
        h_1 = .5*std::pow(x.at(0) - abs_u_k, 2) + 
          (lambda_i/L)*std::log(1.0 + x.at(0)/tuningParameters.theta);
        h_2 = .5*std::pow(x.at(1) - abs_u_k, 2) + 
          (lambda_i/L)*std::log(1.0 + x.at(1)/tuningParameters.theta);
        
        if( h_1 <= h_2 ){
          
          parameters_kp1.at(p) = x.at(0);
          
        }else{
          
          parameters_kp1.at(p) = x.at(1);
          
        }
        
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
      double absPar, lambda_p;
      for(int p = 0; p < parameterValues.n_elem; p ++){
        
        absPar = std::abs(parameterValues.at(p));
        lambda_p = tuningParameters.weights.at(p) * tuningParameters.lambda;
        
        if(absPar <= lambda_p * tuningParameters.theta){
          
          penalty += lambda_p * absPar - std::pow( parameterValues.at(p), 2 ) / (2*tuningParameters.theta);
          
        }else if(absPar > lambda_p * tuningParameters.theta){
          
          penalty += tuningParameters.theta * std::pow(lambda_p, 2) / 2.0;
          
        }else{
          Rcpp::stop("Error while evaluating mcp");
        }
        
      }
      
      return penalty;
    }
    
  };
  
}
#endif
