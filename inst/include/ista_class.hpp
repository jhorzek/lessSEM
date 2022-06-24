#ifndef OPTIMIZERCLASS_H
#define OPTIMIZERCLASS_H
#include <RcppArmadillo.h>
#include "model.hpp"
#include "proximalOperator.hpp"
#include "penalty.hpp"
#include "smoothPenalty.hpp"

// The idea follows ensmallen in that the user supplies a C++ class with
// methods fit and gradients which is used by the optimizer

namespace linr{

struct control{
  const double L0;
  const double eta;
  const int maxIterOut;
  const int maxIterIn;
  const double breakOuter;
  const int verbose;
};

struct fitResults{
  double fit;
  arma::colvec fits;
  bool convergence;
  Rcpp::NumericVector parameterValues;
};

template<typename T> // T is the type of the tuning parameters
inline fitResults ista(model& model_, 
                       Rcpp::NumericVector startingValues,
                       proximalOperator<T>& proximalOperator_, // proximalOperator takes the tuning parameters
                       // as input -> <T>
                       penalty<T>& penalty_, // penalty takes the tuning parameters
                       smoothPenalty<T>& smoothPenalty_, // smoothPenalty takes the tuning parameters
                       // as input -> <T>
                       const T& tuningParameters, // tuning parameters are of type T
                       const control& control_){
  
  Rcpp::StringVector parameterLabels = startingValues.names();
  
  Rcpp::NumericVector parameters_k = Rcpp::clone(startingValues), 
    parameters_kMinus1 = Rcpp::clone(startingValues);
  arma::rowvec parameterChange(startingValues.length());
  arma::rowvec armaGradients(startingValues.length());
  
  double fit_k = model_.fit(startingValues), 
    fit_kMinus1 = model_.fit(startingValues),
    penalty_k = 0.0;
  double penalizedFit_k, penalizedFit_kMinus1;
  Rcpp::NumericVector gradients_k, gradients_kMinus1;
  bool breakInner = false, breakOuter = false;
  arma::mat quadr, parchTimeGrad;
  
  double ridgePenalty = 0.0;
  
  
  penalizedFit_k = fit_k + 
    penalty_.getValue(parameters_k, tuningParameters) + // lasso penalty part
    smoothPenalty_.getValue(parameters_k, tuningParameters); // ridge penalty part
  penalizedFit_kMinus1 = fit_kMinus1 + 
    penalty_.getValue(parameters_k, tuningParameters) + // lasso penalty part
    smoothPenalty_.getValue(parameters_k, tuningParameters); // ridge penalty part
  
  gradients_k = model_.gradients(startingValues) +
    smoothPenalty_.getGradients(parameters_k, tuningParameters); // ridge part
  gradients_kMinus1 = model_.gradients(startingValues) +
    smoothPenalty_.getGradients(parameters_k, tuningParameters); // ridge part
  
  // initialize step size
  double L_kMinus1 = control_.L0, L_k = control_.L0;
  
  arma::colvec fits(control_.maxIterOut);
  fits.fill(arma::datum::nan);
  
  // outer iteration
  double eta = control_.eta;
  for(int outer_iteration = 0; outer_iteration < control_.maxIterOut; outer_iteration ++){
    // TODO: improve initial step size (e.g., using Barzilai Borwein; see GIST)
    if(control_.verbose==-99) Rcpp::Rcout << "Outer iteration " << outer_iteration + 1 << std::endl;
    if(control_.verbose==-99) Rcpp::Rcout << parameters_k << std::endl;
    
    // inner iteration: iterate over all parameters until the convergence criterion is met
    
    for(int inner_iteration = 0; inner_iteration < control_.maxIterIn; inner_iteration ++){
      L_k = std::pow(eta,inner_iteration)*control_.L0;
      if(control_.verbose==-99) Rcpp::Rcout << "L_k : " << L_k << std::endl;
      if(control_.verbose==-99) Rcpp::Rcout << "gradients_k : " << gradients_k << std::endl;
      if(control_.verbose==-99) Rcpp::Rcout << "gradients_kMinus1 : " << gradients_kMinus1 << std::endl;
      
      
      // apply proximal operator to get F(p_L(y))
      
      parameters_k = proximalOperator_.getParameters(
        parameters_kMinus1,
        gradients_kMinus1,
        L_k,
        tuningParameters
      );
      
      parameterChange = parameters_k-parameters_kMinus1;
      if(control_.verbose==-99) 
      {
        Rcpp::Rcout << "parameters_k" << std::endl;
        Rcpp::Rcout << parameters_k << std::endl;
        Rcpp::Rcout << "parameters_kMinus1" << std::endl;
        Rcpp::Rcout << parameters_kMinus1 << std::endl;
        Rcpp::Rcout << "parameterChange" << std::endl;
        Rcpp::Rcout << parameterChange << std::endl;
      }
      
      
      fit_k = model_.fit(parameters_k);
      
      if(control_.verbose==-99)
      {
        Rcpp::Rcout << "ridge gradients: " << smoothPenalty_.getGradients(parameters_k, tuningParameters) << std::endl;
        Rcpp::Rcout << "fit_k : " << fit_k << std::endl;
      }
      
      if(!arma::is_finite(fit_k)) continue;
      
      penalty_k = penalty_.getValue(parameters_k, tuningParameters);
      penalizedFit_k = fit_k + 
        penalty_k +  // lasso part
        smoothPenalty_.getValue(parameters_k, tuningParameters); // ridge part
      
      if(control_.verbose==-99){
        Rcpp::Rcout << "penalizedFit_k : " << penalizedFit_k << std::endl;
        Rcpp::Rcout << "penalizedFit_kMinus1 : " << penalizedFit_kMinus1 << std::endl;
        Rcpp::Rcout << "ridge penalty: " << smoothPenalty_.getValue(parameters_k, tuningParameters) << std::endl;
      }
      
      if(inner_iteration > 0){
        
        quadr = parameterChange*arma::trans(parameterChange);
        armaGradients = Rcpp::as<arma::rowvec>(gradients_kMinus1);
        parchTimeGrad = parameterChange*arma::trans(armaGradients);
        breakInner = penalizedFit_k <= (fit_kMinus1 + parchTimeGrad(0,0) +
          (L_k/2.0)*quadr(0,0) + penalty_k);
        
        // check convergence
        if(control_.verbose==-1){
          Rcpp::Rcout << "penalizedFit_k : " << penalizedFit_k << std::endl;
          Rcpp::Rcout << "fit_kMinus1 : " << fit_kMinus1 << std::endl;
          Rcpp::Rcout << "parchTimeGrad : " << parchTimeGrad(0,0) << std::endl;
          Rcpp::Rcout << "(L_k/2.0)*quadr(0,0) : " << (L_k/2.0)*quadr(0,0) << std::endl;
          Rcpp::Rcout << "penalty_k : " << penalty_k << std::endl;
          Rcpp::Rcout << "penalizedFit_kMinus1 : " << penalizedFit_kMinus1 << std::endl;
        }
        
        if(breakInner) {
          
          // no need to compute gradients if the fit is not improved
          gradients_k = model_.gradients(parameters_k) + 
            smoothPenalty_.getGradients(parameters_k, tuningParameters); // ridge part
          // if any of the gradients is non-finite, we can skip to a 
          // smaller step size
          if(Rcpp::any(!Rcpp::is_finite(gradients_k)))  continue;
          
          if(control_.verbose==-99){
            Rcpp::Rcout << "gradients_k\n: " << gradients_k << std::endl;
            Rcpp::Rcout << "ridge gradients: " << smoothPenalty_.getGradients(parameters_k, tuningParameters) << std::endl;
          }
          
          // if everything worked out fine, we break the inner iteration
          if(control_.verbose==-99) Rcpp::Rcout << "Breaking inner iteration: " << penalizedFit_k << std::endl;
          break;
        }
      }
    }
    
    if(outer_iteration % control_.verbose == 0){
      Rcpp::Rcout << "Fit in iteration outer_iteration " << outer_iteration + 1 << ": " << penalizedFit_k << std::endl;
      Rcpp::Rcout << parameters_k << std::endl;
    } 
    
    if(!breakInner){
      Rcpp::warning("Inner iterations did not improve the fit. Re-running with different eta");
      eta = std::sqrt(control_.eta);
      continue;
    }
    
    // reset eta if it was changed before
    eta = control_.eta;
    
    fits.at(outer_iteration) = penalizedFit_k;
    
    // check outer breaking condition
    if(outer_iteration > 0){
      breakOuter = std::abs(fits.at(outer_iteration) - fits.at(outer_iteration-1)) < control_.breakOuter;
    }
    if(breakOuter) {
      if(control_.verbose==-99) Rcpp::Rcout << "Breaking outer iteration" << std::endl;
      
      break;
    }
    
    penalizedFit_kMinus1 = penalizedFit_k;
    parameters_kMinus1 = parameters_k;
    gradients_kMinus1 = gradients_k;
    L_kMinus1 = L_k;
  }
  
  fitResults fitResults_;
  
  fitResults_.convergence = breakOuter;
  fitResults_.fit = penalizedFit_k;
  fitResults_.fits = fits;
  fitResults_.parameterValues = parameters_k;
  
  return(fitResults_);
}

}// end namespace

#endif