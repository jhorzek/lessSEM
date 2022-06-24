#ifndef OPTIMIZERCLASS_H
#define OPTIMIZERCLASS_H
#include <RcppArmadillo.h>
#include "model.hpp"
#include "proximalOperator.hpp"
#include "penalty.hpp"
#include "smoothPenalty.hpp"

// The design follows ensmallen (https://github.com/mlpack/ensmallen) in that the 
// user supplies a C++ class with methods fit and gradients which is used
// by the optimizer

namespace linr{

struct control{
  const double L0; // L0 controls the step size used in the first iteration
  const double eta; // eta controls by how much the step size changes in the
  // inner iterations with (eta^i)*L, where i is the inner iteration
  const int maxIterOut; // maximal number of outer iterations
  const int maxIterIn; // maximal number of inner iterations
  const double breakOuter; // change in fit required to break the outer iteration
  const int verbose; // if set to a value > 0, the fit every verbose iterations
  // is printed. If set to -99 you will get the debug output which is horribly
  // convoluted
};

struct fitResults{
  double fit; // the final fit value (regularized fit)
  Rcpp::NumericVector fits; // a vector with all fits at the outer iteration
  bool convergence; // was the outer breaking condition met?
  Rcpp::NumericVector parameterValues; // final parameter values
};

// The implementation of ista follows that outlined in 
// Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
// Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
// 183–202. https://doi.org/10.1137/080716542
// see Remark 3.1 on p. 191 (ISTA with backtracking)


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
  if(control_.verbose != 0) Rcpp::Rcout << "Optimizing with ista" << std::endl;
  // we rely heavily on parameter labels because we want to make sure that we 
  // are changing the correct values
  Rcpp::StringVector parameterLabels = startingValues.names();
  
  // prepare parameter vectors
  Rcpp::NumericVector parameters_k = Rcpp::clone(startingValues), 
    parameters_kMinus1 = Rcpp::clone(startingValues);
  // the following elements will be required to judge the breaking condition
  arma::rowvec parameterChange(startingValues.length());
  arma::rowvec armaGradients(startingValues.length());
  arma::mat quadr, parchTimeGrad;
  
  // prepare fit elements
  double fit_k = model_.fit(startingValues), 
    fit_kMinus1 = model_.fit(startingValues),
    penalty_k = 0.0;
  double penalizedFit_k, penalizedFit_kMinus1;
  Rcpp::NumericVector gradients_k, gradients_kMinus1;
  
  double ridgePenalty = 0.0;
  
  penalizedFit_k = fit_k + 
    penalty_.getValue(parameters_k, tuningParameters) + // lasso penalty part
    smoothPenalty_.getValue(parameters_k, tuningParameters); // ridge penalty part
  penalizedFit_kMinus1 = fit_kMinus1 + 
    penalty_.getValue(parameters_kMinus1, tuningParameters) + // lasso penalty part
    smoothPenalty_.getValue(parameters_kMinus1, tuningParameters); // ridge penalty part
  
  // the following vector will save the fits of all iterations:
  Rcpp::NumericVector fits(control_.maxIterOut+1);
  fits.fill(NA_REAL);
  fits.at(0) = penalizedFit_kMinus1;
  
  // prepare gradient elements 
  // NOTE: We combine the gradients of the smooth functions (the log-Likelihood)
  // of the model and the smooth penalty function (e.g., ridge)
  gradients_k = model_.gradients(parameters_k) +
    smoothPenalty_.getGradients(parameters_k, tuningParameters); // ridge part
  gradients_kMinus1 = model_.gradients(parameters_kMinus1) +
    smoothPenalty_.getGradients(parameters_kMinus1, tuningParameters); // ridge part
  
  // breaking flags
  bool breakInner = false, // if true, the inner iteration is exited
    breakOuter = false; // if true, the outer iteration is exited
  
  // initialize step size
  double L_kMinus1 = control_.L0, L_k = control_.L0;
  
  // outer iteration
  for(int outer_iteration = 0; outer_iteration < control_.maxIterOut; outer_iteration ++){    
    if(control_.verbose == -99) Rcpp::Rcout << "Outer iteration " << outer_iteration + 1 << std::endl;
    
    // TODO: improve initial step size (e.g., using Barzilai Borwein; see GIST)
    
    // the gradients will be used by the inner iteration to compute the new 
    // parameters
    gradients_kMinus1 = model_.gradients(parameters_kMinus1) +
      smoothPenalty_.getGradients(parameters_kMinus1, tuningParameters); // ridge part
    armaGradients = Rcpp::as<arma::rowvec>(gradients_kMinus1); // needed in convergence criterion
    
    if(control_.verbose == -99) Rcpp::Rcout << "gradients_kMinus1 : " << gradients_kMinus1 << std::endl;
    
    for(int inner_iteration = 0; inner_iteration < control_.maxIterIn; inner_iteration ++){
      // inner iteration: reduce step size until the convergence criterion is met
      L_k = std::pow(control_.eta, inner_iteration)*L_kMinus1;
      
      if(control_.verbose == -99) Rcpp::Rcout << "L_k : " << L_k << std::endl;
      
      // apply proximal operator to get new parameters for given step size
      parameters_k = proximalOperator_.getParameters(
        parameters_kMinus1,
        gradients_kMinus1,
        L_k,
        tuningParameters
      );
      
      // compute new fit; if this fit is non-finite, we can jump to the next
      // iteration
      fit_k = model_.fit(parameters_k);
      if(control_.verbose == -99)
      {
        Rcpp::Rcout << "fit_k : " << fit_k << std::endl;
      }
      
      if(!arma::is_finite(fit_k)) continue;
      
      // fit_k is only part of the fit we are interested in. We also need
      // the penalty values:
      penalty_k = penalty_.getValue(parameters_k, tuningParameters);
      penalizedFit_k = fit_k + 
        penalty_k +  // lasso part
        smoothPenalty_.getValue(parameters_k, tuningParameters); // ridge part
      
      if(control_.verbose == -99){
        Rcpp::Rcout << "penalizedFit_k : " << penalizedFit_k << std::endl;
        Rcpp::Rcout << "penalizedFit_kMinus1 : " << penalizedFit_kMinus1 << std::endl;
      }
      
      if(!arma::is_finite(penalizedFit_k)) continue;
      
      // to test the convergence criterion, we also have to compute
      // the approximated fit based on the quadratic approximation
      // h(parameters_k) := fit(parameters_k) + 
      // (parameters_k-parameters_kMinus1)*gradients_k^T + 
      // (L/2)*(parameters_k-parameters_kMinus1)^2 + 
      // penalty(parameters_k) 
      parameterChange = parameters_k-parameters_kMinus1;
      quadr = parameterChange*arma::trans(parameterChange); // always positive
      
      parchTimeGrad = parameterChange*arma::trans(armaGradients); // can be 
      // positive or negative
      
      breakInner = penalizedFit_k <= (
        fit_kMinus1 + 
          parchTimeGrad(0,0) +
          (L_k/2.0)*quadr(0,0) + 
          penalty_k
      );
      
      if(control_.verbose == -99){
        Rcpp::Rcout << "penalizedFit_k : " << penalizedFit_k << std::endl;
        Rcpp::Rcout << "fit_kMinus1 : " << fit_kMinus1 << std::endl;
        Rcpp::Rcout << "parchTimeGrad : " << parchTimeGrad << std::endl;
        Rcpp::Rcout << "quadr : " << quadr << std::endl;
        Rcpp::Rcout << "(L_k/2.0)*quadr(0,0) : " << (L_k/2.0)*quadr(0,0) << std::endl;
        Rcpp::Rcout << "penalty_k : " << penalty_k << std::endl;
        Rcpp::Rcout << "penalizedFit_kMinus1 : " << penalizedFit_kMinus1 << std::endl;
      }
      
      if(breakInner) {
        // compute gradients at new position
        gradients_k = model_.gradients(parameters_k) + 
          smoothPenalty_.getGradients(parameters_k, tuningParameters); // ridge part
        
        if(control_.verbose == -99){
          Rcpp::Rcout << "gradients_k\n: " << gradients_k << std::endl;
        }
        
        // if any of the gradients is non-finite, we can skip to a 
        // smaller step size
        if(Rcpp::any(!Rcpp::is_finite(gradients_k)))  continue;
        
        // if everything worked out fine, we break the inner iteration
        if(control_.verbose == -99) Rcpp::Rcout << "Breaking inner iteration" << std::endl;
        break;
        
      }// end break inner
    }// end inner iteration
    
    // print fit info
    if(control_.verbose > 0 && outer_iteration % control_.verbose == 0){
      Rcpp::Rcout << "Fit in iteration outer_iteration " << outer_iteration + 1 << ": " << penalizedFit_k << std::endl;
      Rcpp::Rcout << parameters_k << std::endl;
    } 
    
    if(!breakInner){
      Rcpp::warning("Inner iterations did not improve the fit. Resetting L");
      L_kMinus1 = control_.L0;
      continue;
    }
    
    fits.at(outer_iteration+1) = penalizedFit_k;
    
    // check outer breaking condition
    if(control_.verbose == -99) Rcpp::Rcout << "Comparing " << fits.at(outer_iteration+1) << " to " << fits.at(outer_iteration) << std::endl;
    breakOuter = std::abs(fits.at(outer_iteration+1) - fits.at(outer_iteration)) < control_.breakOuter;
    
    if(breakOuter) {
      if(control_.verbose == -99) Rcpp::Rcout << "Breaking outer iteration" << std::endl;
      break;
    }
    
    // for next iteration: save current values as previous values
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