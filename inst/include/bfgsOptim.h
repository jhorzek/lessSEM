#ifndef BFGSOPTIMCLASS_H
#define BFGSOPTIMCLASS_H
#include <RcppArmadillo.h>
#include "model.h"
#include "fitResults.h"
#include "proximalOperator.h"
#include "ridge.h"
#include "bfgs.h"

// The design follows ensmallen (https://github.com/mlpack/ensmallen) in that the 
// user supplies a C++ class with methods fit and gradients which is used
// by the optimizer.

// Important: This is not the same implementation of BFGS as that used in, for instance,
// optim. Instead, the objective was to create a BFGS optimizer for smooth functions
// which is very similar to the GLMNET optimizer proposed by 
// 1) Friedman, J., Hastie, T., & Tibshirani, R. (2010). 
// Regularization Paths for Generalized Linear Models via Coordinate Descent. 
// Journal of Statistical Software, 33(1), 1–20. https://doi.org/10.18637/jss.v033.i01
// 2) Yuan, G.-X., Chang, K.-W., Hsieh, C.-J., & Lin, C.-J. (2010).
// A Comparison of Optimization Methods and Software for Large-scale 
// L1-regularized Linear Classification. Journal of Machine Learning Research, 11, 3183–3234.
// 3) Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
// An improved GLMNET for l1-regularized logistic regression. 
// The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421


namespace lessSEM{

enum convergenceCriteriaBFGS{
  GLMNET_,
  fitChange_,
  gradients_
};
const std::vector<std::string> convergenceCriteriaBFGS_txt = {
  "GLMNET_",
  "fitChange_",
  "gradients_"
};


struct controlBFGS{
  const arma::mat initialHessian;
  const double stepSize;
  const double sigma;
  const double gamma;
  const int maxIterOut; // maximal number of outer iterations
  const int maxIterIn; // maximal number of inner iterations
  const int maxIterLine;
  const double breakOuter; // change in fit required to break the outer iteration
  const double breakInner;
  const convergenceCriteriaBFGS convergenceCriterion; // this is related to the inner
  // breaking condition. 
  const int verbose; // if set to a value > 0, the fit every verbose iterations
  // is printed.
};

// the line search procedure is the same as used by glmnet
template<typename T> // T is the type of the tuning parameters
inline arma::rowvec bfgsLineSearch(
    model& model_,
    smoothPenalty<T>& smoothPenalty_,
    const arma::rowvec& parameters_kMinus1, 
    const Rcpp::StringVector& parameterLabels,
    const arma::rowvec& direction,
    const double fit_kMinus1, 
    const arma::rowvec& gradients_kMinus1,
    const arma::mat& Hessian_kMinus1,
    
    const T& tuningParameters,
    
    const double stepSize, 
    const double sigma,
    const double gamma, 
    const int maxIterLine,
    const int verbose){
  
  arma::rowvec gradients_k(gradients_kMinus1.n_rows);
  gradients_k.fill(arma::datum::nan);
  arma::rowvec parameters_k(gradients_kMinus1.n_rows);
  parameters_k.fill(arma::datum::nan);
  
  Rcpp::NumericVector randomNumber;
  
  double fit_k; // new fit value of differentiable part
  double p_k; // new penalty value
  double f_k; // new combined fit
  
  // get penalized M2LL for step size 0:
  // Note: we assume that the penalty is already incorporated in the 
  // normal fit function; the following is only in here to demonstrate the 
  // parallels to glmnet
  double pen_0 = 0.0;
  // Note: we see the smooth penalty as part of the smooth
  // objective function and not as part of the non-differentiable
  // penalty
  double f_0 = fit_kMinus1 + pen_0;
  // needed for convergence criterion (see Yuan et al. (2012), Eq. 20)
  // Note: we assume that the penalty is already incorporated in the 
  // normal fit function; the following is only in here to demonstrate the 
  // parallels to glmnet
  double pen_d = 0.0;
  
  double currentStepSize;
  // a step size of >= 1 would result in no change or in an increasing step
  // size
  if(stepSize >= 1) {
    currentStepSize = .9;
  }else{
    currentStepSize = stepSize;
  }
  
  randomNumber = Rcpp::runif(1,0.0,1.0);
  // randomly resetting the step size can help
  // if the optimizer is stuck
  if(randomNumber.at(0) < 0.25){
    Rcpp::NumericVector tmp = Rcpp::runif(1,.5,.99);
    currentStepSize = tmp.at(0);
  }
  
  bool converged = false;
  
  for(int iteration = 0; iteration < maxIterLine; iteration ++){
    
    // set step size
    currentStepSize = std::pow(stepSize, iteration); // starts with 1 and
    // then decreases with each iteration
    
    parameters_k = parameters_kMinus1 + currentStepSize * direction;
    
    fit_k = model_.fit(parameters_k,
                       parameterLabels) + 
                         smoothPenalty_.getValue(parameters_k, 
                                                 parameterLabels, 
                                                 tuningParameters);
    
    if(!arma::is_finite(fit_k)){
      // skip to next iteration and try a smaller step size
      continue;
    }
    
    // compute g(stepSize) = g(x+td) + p(x+td) - g(x) - p(x), 
    // where g is the differentiable part and p the non-differentiable part
    // p(x+td):
    p_k = 0.0;
    
    // g(x+td) + p(x+td)
    f_k = fit_k + p_k;
    
    // test line search criterion. g(stepSize) must show a large enough decrease
    // to be accepted
    // see Equation 20 in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). 
    // An improved GLMNET for l1-regularized logistic regression. 
    // The Journal of Machine Learning Research, 13, 1999–2030. 
    // https://doi.org/10.1145/2020408.2020421
    
    arma::mat compareTo = 
      gradients_kMinus1*arma::trans(direction) + // gradients and direction typically show
      // in the same direction -> positive
      gamma*(direction*Hessian_kMinus1*arma::trans(direction)) + // always positive
      pen_d - pen_0;
    // gamma is set to zero by Yuan et al. (2012)
    // if sigma is 0, no decrease is necessary
    
    converged = f_k - f_0 <=  sigma*currentStepSize*compareTo(0,0);
    
    if(converged){
      // check if gradients can be computed at the new location; 
      // this can often cause issues
      gradients_k = model_.gradients(parameters_k,
                                     parameterLabels);
      
      if(!arma::is_finite(gradients_k)) {
        // go to next iteration and test smaller step size
        continue;
      }
      // else
      break;
    }
    
  }// end line search
  
  if(!converged) Rcpp::warning("Line search did not converge.");
  
  return(parameters_k);
}

template<typename T> // T is the type of the tuning parameters
inline lessSEM::fitResults bfgsOptim(model& model_, 
                               Rcpp::NumericVector startingValuesRcpp,
                               smoothPenalty<T>& smoothPenalty_,
                               const T& tuningParameters, // tuning parameters are of type T
                               const controlBFGS& control_){
  if(control_.verbose != 0) {
    Rcpp::Rcout << "Optimizing with bfgs.\n" <<
      std::endl;
  }
  
  // separate labels and values
  arma::rowvec startingValues = Rcpp::as<arma::rowvec>(startingValuesRcpp);
  const Rcpp::StringVector parameterLabels = startingValuesRcpp.names();
  
  // prepare parameter vectors
  arma::rowvec parameters_k = startingValues, 
    parameters_kMinus1 = startingValues;
  arma::rowvec direction(startingValues.n_elem);
  
  // prepare fit elements
  // fit of the smooth part of the fit function
  double fit_k = model_.fit(parameters_k,
                            parameterLabels) + 
                              smoothPenalty_.getValue(parameters_k,
                                                      parameterLabels,
                                                      tuningParameters);
  double fit_kMinus1 = model_.fit(parameters_kMinus1,
                                  parameterLabels) + 
                                    smoothPenalty_.getValue(parameters_kMinus1,
                                                            parameterLabels,
                                                            tuningParameters);
  // add non-differentiable part -> there is none here
  double penalizedFit_k = fit_k;
  
  double penalizedFit_kMinus1 = fit_kMinus1;
  
  // the following vector will save the fits of all iterations:
  Rcpp::NumericVector fits(control_.maxIterOut+1);
  fits.fill(NA_REAL);
  fits.at(0) = penalizedFit_kMinus1;
  
  // prepare gradient elements 
  // NOTE: We combine the gradients of the smooth functions (the log-Likelihood)
  // of the model and the smooth penalty function (e.g., ridge)
  arma::rowvec gradients_k = model_.gradients(parameters_k,
                                              parameterLabels) +
                                                smoothPenalty_.getGradients(parameters_k, 
                                                                            parameterLabels, 
                                                                            tuningParameters); // ridge part
  arma::rowvec gradients_kMinus1 = model_.gradients(parameters_kMinus1,
                                                    parameterLabels) +
                                                      smoothPenalty_.getGradients(parameters_kMinus1, 
                                                                                  parameterLabels, 
                                                                                  tuningParameters); // ridge part
  
  // prepare Hessian elements
  arma::mat Hessian_k = control_.initialHessian,
    Hessian_kMinus1 = control_.initialHessian;
  
  // breaking flags
  bool breakOuter = false; // if true, the outer iteration is exited
  
  // outer iteration
  for(int outer_iteration = 0; outer_iteration < control_.maxIterOut; outer_iteration ++){

    // check if user wants to stop the computation:
    Rcpp::checkUserInterrupt();
    
    // the gradients will be used by the inner iteration to compute the new 
    // parameters
    gradients_kMinus1 = model_.gradients(parameters_kMinus1, parameterLabels) +
      smoothPenalty_.getGradients(parameters_kMinus1, parameterLabels, tuningParameters); // ridge part
    
    // find step direction -> simple quasi-Newton step
    direction = -arma::trans(arma::solve(Hessian_kMinus1,arma::trans(gradients_kMinus1)));
    
    // find length of step in direction
    parameters_k = bfgsLineSearch(model_, 
                                  smoothPenalty_, 
                                  parameters_kMinus1, 
                                  parameterLabels,
                                  direction,
                                  fit_kMinus1, 
                                  gradients_kMinus1,
                                  Hessian_kMinus1,
                                  
                                  tuningParameters,
                                  
                                  control_.stepSize, 
                                  control_.sigma,
                                  control_.gamma, 
                                  control_.maxIterLine,
                                  control_.verbose);
    
    // get gradients of differentiable part
    gradients_k = model_.gradients(parameters_k,
                                   parameterLabels) +
                                     smoothPenalty_.getGradients(parameters_k,
                                                                 parameterLabels,
                                                                 tuningParameters);
    // fit of the smooth part of the fit function
    fit_k = model_.fit(parameters_k,
                       parameterLabels) + 
                         smoothPenalty_.getValue(parameters_k,
                                                 parameterLabels,
                                                 tuningParameters);
    // add non-differentiable part -> there is none here
    penalizedFit_k = fit_k;
    
    fits.at(outer_iteration+1) = penalizedFit_k;
    
    // print fit info
    if(control_.verbose > 0 && outer_iteration % control_.verbose == 0){
      Rcpp::Rcout << "Fit in iteration outer_iteration " << 
        outer_iteration + 1 << 
          ": " << penalizedFit_k << 
            std::endl;
      Rcpp::Rcout << parameters_k << std::endl;
    } 
    
    // Approximate Hessian using BFGS
    Hessian_k = lessSEM::BFGS(
      parameters_kMinus1,
      gradients_kMinus1,
      Hessian_kMinus1,
      parameters_k,
      gradients_k,
      true,
      .001,
      control_.verbose == -99
    );
    
    // check convergence 
    if(control_.convergenceCriterion == GLMNET_) {
      arma::mat HessDiag = arma::eye(Hessian_k.n_rows,
                                     Hessian_k.n_cols);
      HessDiag.fill(0.0);
      HessDiag.diag() = Hessian_k.diag();
      try{
        breakOuter = max(HessDiag*arma::pow(arma::trans(direction),2)) < control_.breakOuter;
      }catch(...){
        Rcpp::stop("Error while computing convergence criterion");
      }
      
    }
    if(control_.convergenceCriterion == fitChange_){
      try{
        breakOuter = std::abs(fits.at(outer_iteration+1) - 
          fits.at(outer_iteration)) < 
            control_.breakOuter;
      }catch(...){
        Rcpp::stop("Error while computing convergence criterion");
      }
    }
    if(control_.convergenceCriterion == gradients_){
      try{
        arma::rowvec subGradients = gradients_k;
        
        // check if all gradients are below the convergence criterion:
        breakOuter = arma::sum(arma::abs(subGradients) < control_.breakOuter) == 
          subGradients.n_elem;
      }catch(...){
        Rcpp::stop("Error while computing convergence criterion");
      }
    }
    
    if(breakOuter){
      break;
    }
    
    // for next iteration: save current values as previous values
    fit_kMinus1 = fit_k;
    penalizedFit_kMinus1 = penalizedFit_k;
    parameters_kMinus1 = parameters_k;
    gradients_kMinus1 = gradients_k;
    Hessian_kMinus1 = Hessian_k;
    
  }// end outer iteration
  
  if(!breakOuter){
    Rcpp::warning("Outer iterations did not converge");
  }
  
  fitResults fitResults_;
  
  fitResults_.convergence = breakOuter;
  fitResults_.fit = penalizedFit_k;
  fitResults_.fits = fits;
  fitResults_.parameterValues = parameters_k;
  fitResults_.Hessian = Hessian_k;
  
  return(fitResults_);
  
}// end bfgs

}// end namespace

#endif