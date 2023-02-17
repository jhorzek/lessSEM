#ifndef SMOOTHPENALTY_H
#define SMOOTHPENALTY_H
#include <RcppArmadillo.h>

namespace lessSEM{

template<class T>
class smoothPenalty{
public:
  
  virtual double getValue(const arma::rowvec& parameterValues,
                          const Rcpp::StringVector& parameterLabels,
                          const T& tuningParameters);
  virtual arma::rowvec getGradients(const arma::rowvec& parameterValues,
                                    const Rcpp::StringVector& parameterLabels,
                                    const T& tuningParameters);
};

// define some smooth penalties:
template<class T>
class noSmoothPenalty: public smoothPenalty<T>{
public:
  
  double getValue(const arma::rowvec& parameterValues,
                  const Rcpp::StringVector& parameterLabels,
                  const T& tuningParameters) override
                  {
                    return(0.0);
                  };
  arma::rowvec getGradients(const arma::rowvec& parameterValues,
                            const Rcpp::StringVector& parameterLabels,
                            const T& tuningParameters) override
                            {
                              arma::rowvec gradients(parameterValues.n_elem);
                              gradients.fill(0.0);
                              return(gradients);
                            };
};

struct tuningParametersSmoothElasticNet{
  double lambda;
  double alpha;
  double epsilon; // to make function smooth
  arma::rowvec weights;
}; 

class smoothElasticNet: public smoothPenalty<tuningParametersSmoothElasticNet>{
  
public:
  
  double getValue(const arma::rowvec& parameterValues, 
                  const Rcpp::StringVector& parameterLabels,
                  const tuningParametersSmoothElasticNet& tuningParameters) 
  override {
    
    double penalty = 0.0;
    double lambda_i;
    
    for(unsigned int p = 0; p < parameterValues.n_elem; p ++){
      // lasso part:
      lambda_i = tuningParameters.alpha *
        tuningParameters.lambda * 
        tuningParameters.weights.at(p);
      
      penalty += lambda_i * 
        std::sqrt(std::pow(parameterValues.at(p),2) + tuningParameters.epsilon);
      
      // ridge part
      lambda_i = (1.0 - tuningParameters.alpha) *
        tuningParameters.lambda * 
        tuningParameters.weights.at(p);
      penalty += lambda_i * 
        std::pow(parameterValues.at(p),2);
      
    }
    
    return penalty;
  }
  
  arma::rowvec getGradients(const arma::rowvec& parameterValues, 
                            const Rcpp::StringVector& parameterLabels,
                            const tuningParametersSmoothElasticNet& tuningParameters) override {
                              
                              arma::rowvec gradients(parameterValues.n_elem);
                              gradients.fill(0.0);
                              double lambda_i;
                              
                              for(unsigned int p = 0; p < parameterValues.n_elem; p++){
                                
                                // if not regularized: nothing to do here
                                if(tuningParameters.weights.at(p) == 0) continue;
                                
                                // lasso part:
                                lambda_i = tuningParameters.alpha *
                                  tuningParameters.lambda * 
                                  tuningParameters.weights.at(p);
                                
                                gradients.at(p) += lambda_i * parameterValues.at(p) *
                                  (1.0/std::sqrt(std::pow(parameterValues.at(p),2) + 
                                  tuningParameters.epsilon));
                                
                                // ridge part
                                lambda_i = (1.0 - tuningParameters.alpha) *
                                  tuningParameters.lambda * 
                                  tuningParameters.weights.at(p);
                                gradients.at(p) += lambda_i * 
                                  2.0 * 
                                  parameterValues.at(p);
                                
                              }// end for parameter
                              
                              return(gradients);
                              
                            }
  
  
};

}// end namespace
#endif