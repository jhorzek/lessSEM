#ifndef RIDGE_H
#define RIDGE_H
#include <RcppArmadillo.h>
#include "proximalOperator.hpp"
#include "smoothPenalty.hpp"
#include "enet.hpp" // for definition of tuning parameters

namespace linr{

class penaltyRidge: public smoothPenalty<tuningParametersEnet>{
public:
  
  double getValue(const Rcpp::NumericVector& parameterValues, 
                  const tuningParametersEnet& tuningParameters) override 
                  {
                    // if ridge is not used:
                    if(tuningParameters.alpha == 1) return(0.0);
                    
                    // else
                    Rcpp::StringVector parameterLabels = parameterValues.names();
                    
                    double penalty = 0.0;
                    double lambda_i;
                    Rcpp::String parameterLabel;
                    
                    for(int p = 0; p < parameterLabels.length(); p ++){
                      parameterLabel = parameterLabels.at(p);
                      
                      lambda_i = (1.0-tuningParameters.alpha) *
                        tuningParameters.lambda * 
                        tuningParameters.weights[parameterLabel];
                      
                      penalty += lambda_i * std::pow(parameterValues[parameterLabel], 2);
                    }
                    
                    return penalty;
                  }
  
  Rcpp::NumericVector getGradients(const Rcpp::NumericVector& parameterValues, 
                                   const tuningParametersEnet& tuningParameters) override 
                                   {
                                     
                                     Rcpp::StringVector parameterLabels = parameterValues.names();
                                     Rcpp::NumericVector gradients(parameterValues.length());
                                     gradients.names() = parameterLabels;
                                     gradients.fill(0.0);
                                     // if ridge is not used:
                                     if(tuningParameters.alpha == 1) return(gradients);
                                     
                                     // else
                                     
                                     double lambda_i;
                                     Rcpp::String parameterLabel;
                                     
                                     for(int p = 0; p < parameterLabels.length(); p ++){
                                       parameterLabel = parameterLabels.at(p);
                                       
                                       lambda_i = (1.0-tuningParameters.alpha) *
                                         tuningParameters.lambda * 
                                         tuningParameters.weights[parameterLabel];
                                       
                                       gradients[parameterLabel] = lambda_i * 
                                         2 * 
                                         parameterValues[parameterLabel];
                                     }
                                     
                                     return gradients;
                                   } 
  
  
};


}
#endif
