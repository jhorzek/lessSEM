#ifndef ISTACLASS_H
#define ISTACLASS_H
#include <RcppArmadillo.h>
#include "model.h"
#include "proximalOperator.h"
#include "penalty.h"

//#define DEBUG 0

namespace ista{

struct control{
  const int i_k;
  const double L0;
  const double eta;
  const int maxIterOut;
  const int maxIterIn;
  const double breakOuter;
};

struct fitResults{
  double fit;
  bool convergence;
  Rcpp::NumericVector parameterValues;
};

template<typename T>
inline fitResults ista(model& model_, 
                       Rcpp::NumericVector startingValues,
                       proximalOperator<T>& proximalOperator_,
                       penalty<T>& penalty_,
                       const T& tuningParameters,
                       const control& control_){
  
  Rcpp::StringVector parameterLabels = startingValues.names();
  
  Rcpp::NumericVector parameters_k = Rcpp::clone(startingValues), 
    parameters_kMinus1 = Rcpp::clone(startingValues);
  arma::rowvec parameterChange(startingValues.length());
  
  double fit_k = model_.fit(startingValues), 
    fit_kMinus1 = model_.fit(startingValues);
  double penalizedFit_k, penalizedFit_kMinus1;
  Rcpp::NumericVector gradients_k = model_.gradients(startingValues), 
    gradients_kMinus1 = model_.gradients(startingValues);
  bool breakInner = false, breakOuter = false;
  arma::mat quadr;
  
  // initialize step size
  double L = control_.L0*std::pow(control_.eta,0);
  
  arma::rowvec fits(control_.maxIterOut);
  fits.fill(arma::datum::nan);
  
  // outer iteration
  
  for(int outer_iteration = 0; outer_iteration < control_.maxIterOut; outer_iteration ++){
    // TODO: improve initial step size (e.g., using Barzilai Borwein; see GIST)
#ifdef DEBUG
    Rcpp::Rcout << "Outer iteration " << outer_iteration + 1 << std::endl;
    Rcpp::Rcout << parameters_k << std::endl;
#endif
    // inner iteration: iterate over all parameters until the convergence criterion is met
    
    for(int inner_iteration = 0; inner_iteration < control_.maxIterIn; inner_iteration ++){
      L = control_.L0*std::pow(control_.eta,inner_iteration);
      
      // apply proximal operator to get F(p_L(y))
      
      parameters_k = proximalOperator_.getParameters(
        parameters_kMinus1,
        gradients_kMinus1,
        L,
        tuningParameters
      );
      
      parameterChange = parameters_kMinus1 - parameters_k;
#ifdef DEBUG
      Rcpp::Rcout << "parameters_k" << std::endl;
      Rcpp::Rcout << parameters_k << std::endl;
      Rcpp::Rcout << "parameters_kMinus1" << std::endl;
      Rcpp::Rcout << parameters_kMinus1 << std::endl;
      Rcpp::Rcout << "parameterChange" << std::endl;
      Rcpp::Rcout << parameterChange << std::endl;
#endif
      
      fit_k = model_.fit(parameters_k);
      
      if(!arma::is_finite(fit_k)) continue;
      
      penalizedFit_k = fit_k + penalty_.getValue(parameters_k, tuningParameters);
      
      // no need to compute gradients if the fit is not improved
      if(penalizedFit_k >= penalizedFit_kMinus1) continue;
      
      gradients_k = model_.gradients(parameters_k);
      
      // if any of the gradients is non-finite, we can skip to a 
      // smaller step size
      if(Rcpp::any(!Rcpp::is_finite(gradients_k)))  continue;
      
      if(inner_iteration > 0){

        // check convergence
        
        quadr = parameterChange*arma::trans(parameterChange);
        breakInner = penalizedFit_k <= penalizedFit_kMinus1 - 
          (L/2.0)*quadr(0,0);

        if(breakInner) {
#ifdef DEBUG
          Rcpp::Rcout << "Breaking inner iteration: " << fits << std::endl;
#endif
          break;
          }
      }
    }
    fits.at(outer_iteration) = fit_k;
    
    // check outer breaking condition
    if(outer_iteration > 0){
      breakOuter = std::abs(fits.at(outer_iteration) - fits.at(outer_iteration-1)) < control_.breakOuter;
    }
    if(breakOuter) {
#ifdef DEBUG
      Rcpp::Rcout << "Breaking outer iteration: " << fits << std::endl;
#endif
      break;
      }
    
    penalizedFit_kMinus1 = penalizedFit_k;
    parameters_kMinus1 = parameters_k;
    gradients_kMinus1 = gradients_k;
    
  }
  
  fitResults fitResults_;
  
  fitResults_.convergence = breakOuter;
  fitResults_.fit = fit_k;
  fitResults_.parameterValues = parameters_k;
  
  return(fitResults_);
}

}// end namespace

#endif