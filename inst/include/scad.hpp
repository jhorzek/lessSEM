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

// The original software package is by Gong et al. is no longer online, but the files
// are also in the package published by 
// Cheng, W., Zhou, H., & Guo, Z. (2022). A Second Order Algorithm for 
// MCP Regularized Optimization. IEEE Access, 10, 42802–42813. 
// https://doi.org/10.1109/ACCESS.2020.3020834


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
    
    return(lambda * absPar);
    
  }else if((lambda < absPar) & (absPar <= lambda * theta)){
    
    return(
      (-std::pow(par, 2) + 
        2.0 * theta * lambda * absPar - std::pow(lambda,2))/
          (2.0 * (theta - 1.0))
    );
    
  }else if( absPar >  (lambda * theta)){
    
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
    
    double abs_u_k;
    Rcpp::String parameterLabel;
    std::vector<double> x(3, 0.0);
    std::vector<double> h(3, 0.0);
    int sign;
    double thetaXlambda = tuningParameters.theta*tuningParameters.lambda;
    
    for(int p = 0; p < parameterValues.n_elem; p ++)
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
      
      x.at(0) = sign * std::min(
        tuningParameters.lambda, 
        std::max(
          0.0,
          abs_u_k - tuningParameters.lambda/L
        )
      );
      
      x.at(1) = sign * std::min(
        thetaXlambda, 
        std::max(
          tuningParameters.lambda,
          ((L * abs_u_k*(tuningParameters.theta - 1.0)) - 
            thetaXlambda)/
              ( L*(tuningParameters.theta - 2.0) )
        )
      );
      
      x.at(2) = sign * std::max(
        thetaXlambda, 
        abs_u_k
      );
      
      double penalty;
      for(int c = 0; c < 3; c++){
        
        // scad penalty value:
        penalty = scadPenalty(
          x.at(c), 
          tuningParameters.lambda,
          tuningParameters.theta
        );
        
        h.at(c) = .5 * std::pow(x.at(c) - u_k.at(p), 2) + // distance between parameters
          (1.0/L) * penalty; // add penalty value of scad:
        Rcpp::Rcout << "x.at("<<c<<") = " << x.at(c) << std::endl;
        Rcpp::Rcout << "L = " << L << std::endl;
        Rcpp::Rcout << ".5 * std::pow(x.at(c) - u_k.at(p), 2) = "<< .5 * std::pow(x.at(c) - u_k.at(p), 2) << std::endl;
        Rcpp::Rcout << "(1.0/L) * penalty = "<< (1.0/L) * penalty << std::endl;
        
        Rcpp::Rcout << "h.at("<<c<<") = " << h.at(c) << std::endl;
      }
      
      for(int c = 0; c < 3; c ++){
        Rcpp::Rcout << "x.at("<<c<<") = " << x.at(c) << std::endl;
        Rcpp::Rcout << "h.at("<<c<<") = " << h.at(c) << std::endl;
      }

      parameters_kp1.at(p) = x.at(std::distance(std::begin(h), 
                                  std::min_element(h.begin(), h.end())));
      
      Rcpp::Rcout << "Selecting " << std::distance(std::begin(h), 
                                         std::min_element(h.begin(), h.end()))
        << ", x = " << parameters_kp1.at(p) << std::endl;
      
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
    double absPar;
    for(int p = 0; p < parameterValues.n_elem; p ++){
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
