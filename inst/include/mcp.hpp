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

inline double mcpPenalty(const double par, 
                         const double lambda_p,
                         const double theta_p) {
  
  double absPar = std::abs(par);
  
  if(absPar <= (lambda_p * theta_p)){
    return(
      lambda_p * absPar - std::pow( par, 2 ) / (2.0 * theta_p)
    );
    
  }else if(absPar > (lambda_p * theta_p)){
    
    return(
      theta_p * std::pow(lambda_p, 2) / 2.0
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
    
    double lambda_p, theta_p, abs_u_k;
    Rcpp::String parameterLabel;
    std::vector<double> C(3, 0.0);
    std::vector<double> f(3, 0.0);
    std::vector<double> x(2, 0.0);
    std::vector<double> h(2, 0.0);
    int sign, nval;
    
    for(int p = 0; p < parameterValues.n_elem; p ++)
    {
#ifdef printit
      Rcpp::Rcout << "####### u_k for " << parameterLabels.at(p) << " = " <<  u_k.at(p) << "#######" << std::endl;
#endif
      lambda_p = tuningParameters.weights.at(p) * tuningParameters.lambda;
      theta_p = tuningParameters.weights.at(p) * tuningParameters.theta;
      
      sign = (u_k.at(p) > 0);
      if(u_k.at(p) < 0) sign = -1;
      
      abs_u_k = std::abs(u_k.at(p));
      
      // compute z first; needed for x_1
      C.at(1) = lambda_p * theta_p;
      
      if(theta_p - 1 != 0){
        
        C.at(2) = std::min(
          C.at(1),
          std::max(
            0.0,
            (theta_p *(L*abs_u_k - lambda_p)) / 
              (L * (theta_p - 1.0))
          )
        );
        
        nval = 3;
        
      }else{
        
        // only the first two values are valid otherwise
        nval = 2;
        
      }
      
#ifdef printit
      Rcpp::Rcout << " nval = " << nval << "\n";
#endif
      for(int w = 0; w < nval; w++){
        f.at(w) = .5*std::pow( C.at(w) - abs_u_k , 2) + 
          (lambda_p/L) * C.at(w) - 
          std::pow(C.at(w),2)/(2*theta_p);
#ifdef printit
        Rcpp::Rcout << " w = " << w << ", C.at(w) = " << C.at(w)  << ", f = " << f.at(w) << "\n";
#endif
      }
#ifdef printit 
      Rcpp::Rcout << " min is achieved with C " << C.at(std::distance(std::begin(f), 
                                              std::min_element(f.begin(), f.end() - (3-nval)))) << "\n";
#endif
      x.at(0) = sign * C.at(std::distance(std::begin(f), 
                            std::min_element(f.begin(), f.end() - (3-nval))));
      
      x.at(1) = sign * std::max( theta_p * lambda_p, abs_u_k );
#ifdef printit 
      Rcpp::Rcout << "x_1 = " <<  x.at(0) << ", x_2 = " <<  x.at(1) << std::endl;
#endif     
      double absPar, penalty;
      for(int c = 0; c < 2; c++){
        
        // mcp penalty value:
        penalty = mcpPenalty(x.at(c), 
                              lambda_p,
                              theta_p);
        
#ifdef printit 
        Rcpp::Rcout << "penalty = " <<  penalty << std::endl;
#endif 
        h.at(c) = .5*std::pow( x.at(c) - u_k.at(p) , 2) + 
          (1.0/L) * penalty;
#ifdef printit 
        Rcpp::Rcout << " c = " << c << ", h = " << h.at(c) << "\n";   
#endif 
      }
      
      if( h.at(0) <= h.at(1) ){
#ifdef printit 
        Rcpp::Rcout << " h.at(0) -> " << x.at(0) << "\n";
#endif 
        parameters_kp1.at(p) = x.at(0);
        
      }else{
#ifdef printit 
        Rcpp::Rcout << " h.at(1) -> " << x.at(1) << "\n";
#endif 
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
    double absPar, lambda_p, theta_p;
    for(int p = 0; p < parameterValues.n_elem; p ++){
      
      lambda_p = tuningParameters.weights.at(p) * tuningParameters.lambda;
      theta_p = tuningParameters.weights.at(p) * tuningParameters.theta;
      
      // mcp penalty value:
      penalty += mcpPenalty(parameterValues.at(p), 
                             lambda_p,
                             theta_p);
      
    }
    
    return penalty;
  }
  
};

}
#endif
