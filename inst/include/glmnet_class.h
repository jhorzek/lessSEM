#ifndef GLMNETCLASS_H
#define GLMNETCLASS_H
#include <RcppArmadillo.h>
#include "model.h"
#include "fitResults.h"
#include "lasso.h"
#include "ridge.h"
#include "enet.h"
#include "bfgs.h"

// The design follows ensmallen (https://github.com/mlpack/ensmallen) in that the 
// user supplies a C++ class with methods fit and gradients which is used
// by the optimizer

// The implementation of GLMNET follows that outlined in 
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

enum convergenceCriteriaGlmnet{
  GLMNET,
  fitChange,
  gradients
};
const std::vector<std::string> convergenceCriteriaGlmnet_txt = {
  "GLMNET",
  "fitChange",
  "gradients"
};


struct controlGLMNET{
  arma::mat initialHessian;
  double stepSize;
  double sigma;
  double gamma;
  int maxIterOut; // maximal number of outer iterations
  int maxIterIn; // maximal number of inner iterations
  int maxIterLine;
  double breakOuter; // change in fit required to break the outer iteration
  double breakInner;
  convergenceCriteriaGlmnet convergenceCriterion; // this is related to the inner
  // breaking condition. 
  int verbose; // if set to a value > 0, the fit every verbose iterations
  // is printed.
};

inline controlGLMNET controlGlmnetDefault(){
  arma::mat initialHessian(1,1);
  initialHessian.fill(1.0);
  controlGLMNET defaultIs ={
    initialHessian, // initial hessian will be diagonal with 100 as diagonal values
    .9, // stepSize;
    1e-5, // sigma
    0, // gamma
    1000, //maxIterOut; // maximal number of outer iterations
    1000, // maxIterIn; // maximal number of inner iterations
    500, // maxIterLine;
    1e-8, // breakOuter; // change in fit required to break the outer iteration
    1e-10, // breakInner;
    fitChange, // convergenceCriterion; // this is related to the inner
    // breaking condition. 
    0 // verbose; // if set to a value > 0, the fit every verbose iterations
    // is printed. 
  };
  return(defaultIs);
}

inline arma::rowvec glmnetInner(const arma::rowvec& parameters_kMinus1,
                                const arma::rowvec& gradients_kMinus1,
                                const arma::mat& Hessian,
                                const arma::rowvec& lambda,
                                const arma::rowvec& alpha,
                                const arma::rowvec& weights,
                                const int maxIterIn,
                                const double breakInner,
                                const int verbose
){
  arma::rowvec stepDirection = parameters_kMinus1;
  stepDirection.fill(0.0);
  arma::rowvec z = parameters_kMinus1;
  z.fill(0.0);
  arma::rowvec parameters_k = parameters_kMinus1;
  parameters_k.fill(0.0);
  arma::colvec hessianXdirection, 
  HessTimesZ(Hessian.n_rows, arma::fill::zeros);
  arma::mat HessDiag(Hessian.n_rows, Hessian.n_cols, arma::fill::zeros),
  dp_k(1,1, arma::fill::zeros), 
  d2p_k(1,1, arma::fill::zeros), 
  z_j(1,1, arma::fill::zeros),
  newParameter(1,1, arma::fill::zeros), 
  zChange(1,1, arma::fill::zeros);
  
  HessDiag.diag() = Hessian.diag();

  // the order in which parameters are updated should be random
  Rcpp::NumericVector randOrder(stepDirection.n_elem);
  Rcpp::NumericVector sampleFrom(stepDirection.n_elem);
  for(unsigned int i = 0; i < stepDirection.n_elem; i++) sampleFrom.at(i) = i;
  
  for(int it = 0; it < maxIterIn; it++){
    
    // reset direction z
    z.fill(arma::fill::zeros); 
    //z_old.fill(arma::fill::zeros);
    
    // iterate over parameters in random order
    randOrder = Rcpp::sample(sampleFrom, stepDirection.n_elem,false);
    
    for(unsigned int p = 0; p < stepDirection.n_elem; p++){
      
      // compute derivative elements:
      hessianXdirection = Hessian*arma::trans(stepDirection);
      dp_k = gradients_kMinus1.col(randOrder.at(p)) + hessianXdirection.row(randOrder.at(p));
      d2p_k = Hessian.row(randOrder.at(p)).col(randOrder.at(p));
      
      // if the parameter is regularized:
      if(weights.at(randOrder.at(p)) != 0){
        newParameter = parameters_kMinus1.at(randOrder.at(p));
        
        // adjust stepDirection for regularized parameters
        if((dp_k(0,0)-alpha.col(randOrder.at(p))(0,0)*lambda.col(randOrder.at(p))(0,0)*weights.col(randOrder.at(p))(0,0)) >= 
           (d2p_k(0,0)*(newParameter(0,0) + stepDirection.col(randOrder.at(p))(0,0)))){
          
          // condition 1
          z_j = -(dp_k-alpha.col(randOrder.at(p))(0,0)*lambda.col(randOrder.at(p))(0,0)*weights.col(randOrder.at(p)))/(d2p_k);
          z.col(randOrder.at(p)) = z_j(0,0);
          stepDirection.col(randOrder.at(p)) += z_j(0,0);
          
        }else if((dp_k(0,0)+alpha.col(randOrder.at(p))(0,0)*lambda.col(randOrder.at(p))(0,0)*weights.col(randOrder.at(p))(0,0)) <=
          (d2p_k(0,0)*(newParameter(0,0) + stepDirection.col(randOrder.at(p))(0,0)))){
          
          // condition 2
          z_j = -(dp_k+alpha.col(randOrder.at(p))(0,0)*lambda.col(randOrder.at(p))(0,0)*weights.col(randOrder.at(p)))/(d2p_k);
          z.col(randOrder.at(p)) = z_j(0,0);
          stepDirection.col(randOrder.at(p)) = stepDirection.col(randOrder.at(p)) + z_j(0,0);
          
        }else{
          
          // condition 3
          z_j = -(newParameter+stepDirection.col(randOrder.at(p)));
          z.col(randOrder.at(p)) = z_j(0,0);
          stepDirection.col(randOrder.at(p)) = stepDirection.col(randOrder.at(p))+z_j(0,0);
          
        }
      }else{
        
        // if not regularized: coordinate descent with newton direction
        z_j = -dp_k/d2p_k;
        z.col(randOrder.at(p)) = z_j(0,0);
        stepDirection.col(randOrder.at(p)) = stepDirection.col(randOrder.at(p))+z_j(0,0);
        
      }
    }
    
    // check inner stopping criterion:
    HessTimesZ = HessDiag*arma::pow(arma::trans(z),2);
    
    // zChange = z*arma::trans(z_old);
    // if(useMultipleConvergenceCriteria & (zChange(0,0) < breakInner)){
    //   break;
    // }
    
    if(HessTimesZ.max() < breakInner){
      break;
    }
    // z_old = z;
  }
  
  return(stepDirection);
  
}


inline arma::rowvec glmnetLineSearch(
    model& model_,
    penaltyLASSOGlmnet& penalty_,
    penaltyRidgeGlmnet& smoothPenalty_, 
    const arma::rowvec& parameters_kMinus1, 
    const Rcpp::StringVector& parameterLabels,
    const arma::rowvec& direction,
    const double fit_kMinus1, 
    const arma::rowvec& gradients_kMinus1,
    const arma::mat& Hessian_kMinus1,
    
    const tuningParametersEnetGlmnet& tuningParameters,
    
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
  
  double pen_0 = penalty_.getValue(parameters_kMinus1, 
                                   parameterLabels,
                                   tuningParameters);
  // Note: we see the smooth penalty as part of the smooth
  // objective function and not as part of the non-differentiable
  // penalty
  double f_0 = fit_kMinus1 + pen_0;
  // needed for convergence criterion (see Yuan et al. (2012), Eq. 20)
  double pen_d = penalty_.getValue(parameters_kMinus1 + direction, 
                                   parameterLabels,
                                   tuningParameters);
  
  double currentStepSize;
  // a step size of >= 1 would result in no change or in an increasing step
  // size
  if(stepSize >= 1) {
    currentStepSize = .9;
  }else{
    currentStepSize = stepSize;
  }
  
  randomNumber = Rcpp::runif(1,0.0,1.0);
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
    p_k = penalty_.getValue(parameters_k, 
                            parameterLabels, 
                            tuningParameters);
    
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
  
  return(parameters_k);
}


inline lessSEM::fitResults glmnet(model& model_, 
                                  Rcpp::NumericVector startingValuesRcpp,
                                  penaltyLASSOGlmnet& penalty_,
                                  penaltyRidgeGlmnet& smoothPenalty_, 
                                  const tuningParametersEnetGlmnet tuningParameters,
                                  const controlGLMNET& control_ = controlGlmnetDefault())
{
  
  if(control_.verbose != 0) {
    Rcpp::Rcout << "Optimizing with glmnet.\n" <<
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
  // add non-differentiable part
  double penalizedFit_k = fit_k +
    penalty_.getValue(parameters_k,
                      parameterLabels,
                      tuningParameters);
  
  double penalizedFit_kMinus1 = fit_kMinus1 +
    penalty_.getValue(parameters_kMinus1,
                      parameterLabels,
                      tuningParameters);
  
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
  arma::mat Hessian_k(startingValues.n_elem,startingValues.n_elem, arma::fill::zeros), 
  Hessian_kMinus1(startingValues.n_elem,startingValues.n_elem, arma::fill::zeros);
  if((control_.initialHessian.n_cols == 1) && (control_.initialHessian.n_rows == 1)){
    // Hessian comes from default initializer and has to be redefined
    double hessianValue = control_.initialHessian(0,0);
    Hessian_k.diag().fill(hessianValue);
    Hessian_kMinus1.diag().fill(hessianValue);
    
  }else{
    Hessian_k = control_.initialHessian;
    Hessian_kMinus1 = control_.initialHessian;
  }
  
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
    
    // find step direction
    direction = glmnetInner(parameters_kMinus1,
                            gradients_kMinus1,
                            Hessian_kMinus1,
                            tuningParameters.lambda,
                            tuningParameters.alpha,
                            tuningParameters.weights,
                            control_.maxIterIn,
                            control_.breakInner,
                            control_.verbose);
    
    // find length of step in direction
    parameters_k = glmnetLineSearch(model_, 
                                    penalty_,
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
    // add non-differentiable part
    penalizedFit_k = fit_k +
      penalty_.getValue(parameters_k,
                        parameterLabels,
                        tuningParameters);
    
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
    if(control_.convergenceCriterion == GLMNET) {
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
    if(control_.convergenceCriterion == fitChange){
      try{
        breakOuter = std::abs(fits.at(outer_iteration+1) - 
          fits.at(outer_iteration)) < 
            control_.breakOuter;
      }catch(...){
        Rcpp::stop("Error while computing convergence criterion");
      }
    }
    if(control_.convergenceCriterion == gradients){
      try{
        arma::rowvec subGradients = penalty_.getSubgradients(
          parameters_k,
          gradients_k,
          tuningParameters
        );
        
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
  
}// end glmnet

}// end namespace

#endif
