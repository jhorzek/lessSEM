#ifndef ISTACLASS_H
#define ISTACLASS_H
#include <RcppArmadillo.h>
#include "model.hpp"
#include "fitResults.hpp"
#include "proximalOperator.hpp"
#include "penalty.hpp"
#include "smoothPenalty.hpp"

// The design follows ensmallen (https://github.com/mlpack/ensmallen) in that the 
// user supplies a C++ class with methods fit and gradients which is used
// by the optimizer

// The implementation of ista follows that outlined in 
// Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding 
// Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
// 183–202. https://doi.org/10.1137/080716542
// see Remark 3.1 on p. 191 (ISTA with backtracking)

// GIST can be found in 
// Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). 
// A General Iterative Shrinkage and Thresholding Algorithm for Non-convex 
// Regularized Optimization Problems. Proceedings of the 30th International 
// Conference on Machine Learning, 28(2)(2), 37–45.

namespace linr{

enum convCritInnerIsta{
  istaCrit,
  gistCrit
};
const std::vector<std::string> convCritInnerIsta_txt = {
  "istaCrit",
  "gistCrit"
};

enum stepSizeInheritance{
  initial,
  istaStepInheritance,
  barzilaiBorwein,
  stochasticBarzilaiBorwein
};

const std::vector<std::string> stepSizeInheritance_txt = {
  "initial",
  "istaStepInheritance",
  "barzilaiBorwein",
  "stochasticBarzilaiBorwein"
};

struct control{
  const double L0; // L0 controls the step size used in the first iteration
  const double eta; // eta controls by how much the step size changes in the
  // inner iterations with (eta^i)*L, where i is the inner iteration
  const int maxIterOut; // maximal number of outer iterations
  const int maxIterIn; // maximal number of inner iterations
  const double breakOuter; // change in fit required to break the outer iteration
  const convCritInnerIsta convCritInner; // this is related to the inner
  // breaking condition. 
  // ista, as presented by Beck & Teboulle (2009); see Remark 3.1 on p. 191 (ISTA with backtracking)
  // gist, as presented by Gong et al. (2013) (Equation 3)
  // 
  const double sigma; // sigma in (0,1) is used by the gist convergence criterion. larger
  // simga enforce larger improvement in fit
  const stepSizeInheritance stepSizeIn; // how should step sizes be carried forward?
  // from iteration to iteration? 
  // initial resets the step size to L0 in each iteration
  // istaStepInheritance takes the previous step size as initial value for the 
  // next iteration
  // barzilaiBorwein uses the Barzilai-Borwein procedure
  // stochasticBarzilaiBorwein uses the Barzilai-Borwein procedure, but sometimes
  // resets the step size; this can help when the optimizer is caught in a bad spot.
  const int verbose; // if set to a value > 0, the fit every verbose iterations
  // is printed. If set to -99 you will get the debug output which is horribly
  // convoluted
};

template<typename T> // T is the type of the tuning parameters
inline linr::fitResults ista(model& model_, 
                       Rcpp::NumericVector startingValuesRcpp,
                       proximalOperator<T>& proximalOperator_, // proximalOperator takes the tuning parameters
                       // as input -> <T>
                       penalty<T>& penalty_, // penalty takes the tuning parameters
                       smoothPenalty<T>& smoothPenalty_, // smoothPenalty takes the tuning parameters
                       // as input -> <T>
                       const T& tuningParameters, // tuning parameters are of type T
                       const control& control_){
  if(control_.verbose != 0) {
    Rcpp::Rcout << "Optimizing with ista.\n" <<
      "Using " << convCritInnerIsta_txt.at(control_.convCritInner) << " as inner convergence criterion\n" <<
        "Using " << stepSizeInheritance_txt.at(control_.stepSizeIn) << " as step size inheritance\n" << 
          std::endl;
  }
  // separate labels and values
  const arma::rowvec startingValues = Rcpp::as<arma::rowvec>(startingValuesRcpp);
  const Rcpp::StringVector parameterLabels = startingValuesRcpp.names();
  
  // prepare parameter vectors
  arma::rowvec parameters_k = startingValues, 
    parameters_kMinus1 = startingValues;
  // the following elements will be required to judge the breaking condition
  arma::rowvec parameterChange(startingValues.n_elem);
  arma::rowvec gradientChange(startingValues.n_elem); // necessary for Barzilai Borwein
  arma::mat quadr, parchTimeGrad;
  
  // prepare fit elements
  double fit_k = model_.fit(startingValues, parameterLabels) +
    smoothPenalty_.getValue(parameters_k, parameterLabels, tuningParameters), // ridge penalty part
    fit_kMinus1 = model_.fit(startingValues, parameterLabels) +
      smoothPenalty_.getValue(parameters_kMinus1, parameterLabels, tuningParameters), // ridge penalty part,
    penalty_k = 0.0;
  double penalizedFit_k, penalizedFit_kMinus1;
  arma::rowvec gradients_k, gradients_kMinus1;
  
  double ridgePenalty = 0.0;
  
  penalizedFit_k = fit_k + 
    penalty_.getValue(parameters_k, parameterLabels, tuningParameters); // lasso penalty part
    
  penalizedFit_kMinus1 = fit_kMinus1 + 
    penalty_.getValue(parameters_kMinus1, parameterLabels, tuningParameters); // lasso penalty part
  
  // the following vector will save the fits of all iterations:
  Rcpp::NumericVector fits(control_.maxIterOut+1);
  fits.fill(NA_REAL);
  fits.at(0) = penalizedFit_kMinus1;
  
  // prepare gradient elements 
  // NOTE: We combine the gradients of the smooth functions (the log-Likelihood)
  // of the model and the smooth penalty function (e.g., ridge)
  gradients_k = model_.gradients(parameters_k, parameterLabels) +
    smoothPenalty_.getGradients(parameters_k, parameterLabels, tuningParameters); // ridge part
  gradients_kMinus1 = model_.gradients(parameters_kMinus1, parameterLabels) +
    smoothPenalty_.getGradients(parameters_kMinus1, parameterLabels, tuningParameters); // ridge part
  
  // breaking flags
  bool breakInner = false, // if true, the inner iteration is exited
    breakOuter = false; // if true, the outer iteration is exited
  
  // initialize step size
  double L_kMinus1 = control_.L0, L_k = control_.L0;
  
  // outer iteration
  for(int outer_iteration = 0; outer_iteration < control_.maxIterOut; outer_iteration ++){    
    if(control_.verbose == -99) Rcpp::Rcout << "Outer iteration " << outer_iteration + 1 << std::endl;
    
    for(int inner_iteration = 0; inner_iteration < control_.maxIterIn; inner_iteration ++){
      // inner iteration: reduce step size until the convergence criterion is met
      L_k = std::pow(control_.eta, inner_iteration)*L_kMinus1;
      
      if(control_.verbose == -99) Rcpp::Rcout << "L_k : " << L_k << std::endl;
      
      // apply proximal operator to get new parameters for given step size
      parameters_k = proximalOperator_.getParameters(
        parameters_kMinus1,
        gradients_kMinus1,
        parameterLabels,
        L_k,
        tuningParameters
      );
      
      // compute new fit; if this fit is non-finite, we can jump to the next
      // iteration
      fit_k = model_.fit(parameters_k, parameterLabels) +
        smoothPenalty_.getValue(parameters_k, parameterLabels, tuningParameters); // ridge penalty part
      if(control_.verbose == -99)
      {
        Rcpp::Rcout << "fit_k : " << fit_k << std::endl;
      }
      
      if(!arma::is_finite(fit_k)) continue;
      
      // fit_k is only part of the fit we are interested in. We also need
      // the penalty values:
      penalty_k = penalty_.getValue(parameters_k, 
                                    parameterLabels, 
                                    tuningParameters);
      
      penalizedFit_k = fit_k + 
        penalty_k;  // lasso part
      
      if(control_.verbose == -99){
        Rcpp::Rcout << "penalizedFit_k : " << penalizedFit_k << std::endl;
        Rcpp::Rcout << "penalizedFit_kMinus1 : " << penalizedFit_kMinus1 << std::endl;
      }
      
      if(!arma::is_finite(penalizedFit_k)) continue;
      
      // to test the convergence criterion, we offer different criteria
      
      if(control_.convCritInner == istaCrit) {
        // ISTA:
        // The approximated fit based on the quadratic approximation
        // h(parameters_k) := fit(parameters_k) + 
        // (parameters_k-parameters_kMinus1)*gradients_k^T + 
        // (L/2)*(parameters_k-parameters_kMinus1)^2 + 
        // penalty(parameters_k) 
        // is compared to the exact fit
        parameterChange = parameters_k-parameters_kMinus1;
        quadr = parameterChange*arma::trans(parameterChange); // always positive
        parchTimeGrad = parameterChange*arma::trans(gradients_kMinus1); // can be 
        // positive or negative
        
        if(control_.verbose == -99){
          Rcpp::Rcout << "penalizedFit_k : " << penalizedFit_k  << " vs " << fit_kMinus1 + 
            parchTimeGrad(0,0) +
            (L_k/2.0)*quadr(0,0) + 
            penalty_k << std::endl;
        }
        
        breakInner = penalizedFit_k <= (
          fit_kMinus1 + 
            parchTimeGrad(0,0) +
            (L_k/2.0)*quadr(0,0) + 
            penalty_k
        );
        
      }else if(control_.convCritInner == gistCrit){
        
        // GIST:
        // the exact fit is compared to 
        // h(parameters_k) := fit(parameters_k) + 
        // penalty(parameters_kMinus1) +
        // L*(sigma/2)*(parameters_k-parameters_kMinus1)^2 + 
        // 
        parameterChange = parameters_k-parameters_kMinus1;
        quadr = parameterChange*arma::trans(parameterChange); // always positive
        
        if(control_.verbose == -99){
          Rcpp::Rcout << "penalizedFit_k : " << penalizedFit_k  << " vs " << penalizedFit_kMinus1 -
            L_k*(control_.sigma/2.0)*quadr(0,0) << std::endl;
        }
        
        breakInner = penalizedFit_k <= (
          penalizedFit_kMinus1 -
            L_k*(control_.sigma/2.0)*quadr(0,0)
        );
      }
      
      if(breakInner) {
        // compute gradients at new position
        gradients_k = model_.gradients(parameters_k, 
                                       parameterLabels) + 
          smoothPenalty_.getGradients(parameters_k, 
                                      parameterLabels, 
                                      tuningParameters); // ridge part
        
        if(control_.verbose == -99){
          Rcpp::Rcout << "gradients_k\n: " << gradients_k << std::endl;
        }
        
        // if any of the gradients is non-finite, we can skip to a 
        // smaller step size
        if(!arma::is_finite(gradients_k))  continue;
        
        // if everything worked out fine, we break the inner iteration
        if(control_.verbose == -99) Rcpp::Rcout << "Breaking inner iteration" << std::endl;
        break;
        
      }// end break inner
    }// end inner iteration
    
    // print fit info
    if(control_.verbose > 0 && outer_iteration % control_.verbose == 0){
      Rcpp::Rcout << "Fit in iteration outer_iteration " << 
        outer_iteration + 1 << 
          ": " << penalizedFit_k << 
            std::endl;
      Rcpp::Rcout << parameters_k << std::endl;
    } 
    
    if(!breakInner){
      Rcpp::warning("Inner iterations did not improve the fit --> resetting L.");
      L_kMinus1 = control_.L0;
      continue;
    }
    
    gradients_k = model_.gradients(parameters_k, 
                                   parameterLabels) + 
                                     smoothPenalty_.getGradients(parameters_k, 
                                                                 parameterLabels, 
                                                                 tuningParameters); // ridge part
    
    fits.at(outer_iteration+1) = penalizedFit_k;
    
    // check outer breaking condition
    if(control_.verbose == -99) Rcpp::Rcout << "Comparing " << fits.at(outer_iteration+1) << " to " << fits.at(outer_iteration) << std::endl;
    breakOuter = std::abs(fits.at(outer_iteration+1) - fits.at(outer_iteration)) < control_.breakOuter;
    
    if(breakOuter) {
      if(control_.verbose == -99) Rcpp::Rcout << "Breaking outer iteration" << std::endl;
      break;
    }
    
    // define new initial step size
    if(control_.stepSizeIn == initial){
      
      L_kMinus1 = control_.L0;
      
    }else if(control_.stepSizeIn == barzilaiBorwein || 
      control_.stepSizeIn == stochasticBarzilaiBorwein){
      
      parameterChange = parameters_k - parameters_kMinus1;
      gradientChange = gradients_k - gradients_kMinus1;

      quadr = (parameterChange)*arma::trans(parameterChange);
      parchTimeGrad = parameterChange*arma::trans(gradientChange);
      if(control_.verbose == -99) Rcpp::Rcout << "quadr: " << quadr << std::endl;
      if(control_.verbose == -99) Rcpp::Rcout << "parchTimeGrad: " << parchTimeGrad << std::endl;
      L_kMinus1 = parchTimeGrad(0,0)/quadr(0,0);
      if(control_.verbose == -99) Rcpp::Rcout << "L_kMinus1 after BB: " << L_kMinus1 << std::endl;
      if(L_kMinus1 < 1e-10 || L_kMinus1 > 1e10) L_kMinus1 = control_.L0;
      
      if(control_.stepSizeIn == stochasticBarzilaiBorwein && 
         rand() % 100 < 25) {
        if(control_.verbose == -99) Rcpp::Rcout << "Resetting L_kMinus1 randomly " << std::endl;
        L_kMinus1 = control_.L0;// reset with 25% probability
        }
      
    }else if(control_.stepSizeIn == istaStepInheritance){
      
      L_kMinus1 = L_k;
      
    }else{
      Rcpp::stop("Unknown step inheritance.");
    }
    
    // for next iteration: save current values as previous values
    fit_kMinus1 = fit_k;
    penalizedFit_kMinus1 = penalizedFit_k;
    parameters_kMinus1 = parameters_k;
    gradients_kMinus1 = gradients_k;
    
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