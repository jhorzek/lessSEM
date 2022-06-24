#ifndef GLMNETCLASS_H
#define GLMNETCLASS_H
#include <RcppArmadillo.h>
#include "model.hpp"
#include "proximalOperator.hpp"
#include "lasso.hpp"
#include "ridge.hpp"
#include "enet.hpp"

// The design follows ensmallen (https://github.com/mlpack/ensmallen) in that the 
// user supplies a C++ class with methods fit and gradients which is used
// by the optimizer

namespace linr{

enum convCritInnerGlmnet{
  glmnetCrit,
  parameterChange
};
const std::vector<std::string> convCritInnerGlmnet_txt = {
  "glmnetCrit",
  "parameterChange"
};

struct controlGLMNET{
  const arma::mat initialHessian;
  const double stepSize;
  const double sigma;
  const double gamma;
  const int maxIterOut; // maximal number of outer iterations
  const int maxIterIn; // maximal number of inner iterations
  const int maxIterLine;
  const double breakOuter; // change in fit required to break the outer iteration
  const double breakInner;
  const convCritInnerGlmnet convCritInner; // this is related to the inner
  // breaking condition. 
  const int verbose; // if set to a value > 0, the fit every verbose iterations
  // is printed. If set to -99 you will get the debug output which is horribly
  // convoluted
};

// The implementation of glmnet follows that outlined in 
// Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1). https://doi.org/10.18637/jss.v033.i01
// and Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421

inline arma::colvec glmnetInner(const arma::rowvec& parameters_kMinus1,
                                const arma::colvec& gradients_kMinus1,
                                const arma::mat& Hessian,
                                const double lambda,
                                const double alpha,
                                const arma::colvec& weights,
                                const int maxIterIn,
                                const double breakInner
){
  arma::colvec stepDirection = Rcpp::clone(parameters_kMinus1);
  stepDirection.fill(0.0);
  arma::colvec z = Rcpp::clone(parameters_kMinus1);
  z.fill(0.0);
  arma::colvec parameters_k = Rcpp::clone(parameters_kMinus1);
  parameters_k.fill(0.0);
  arma::colvec hessianXdirection, 
  HessTimesZ(Hessian.n_rows, arma::fill::zeros);
  arma::mat HessDiag(Hessian.n_rows, Hessian.n_cols, arma::fill::zeros),
  dp_k(1,1, arma::fill::zeros), 
  d2p_k(1,1, arma::fill::zeros), 
  z_j(1,1, arma::fill::zeros),
  newParameter(1,1, arma::fill::zeros), 
  zChange(1,1, arma::fill::zeros);
  
  std::vector<int> randOrder(stepDirection.n_elem);
  for(int i = 0 ; i < stepDirection.n_elem; i++) randOrder.at(i) = i;
  
  for(int it = 0; it < maxIterIn; it++){
    
    // reset direction z
    z.fill(arma::fill::zeros); 
    //z_old.fill(arma::fill::zeros);
    // iterate over parameters in random order
    std::random_shuffle(randOrder.begin(), randOrder.end());
    
    for(int p = 0; p < stepDirection.n_elem; p++){
      
      // compute derivative elements:
      hessianXdirection = Hessian*stepDirection;
      dp_k = gradients_kMinus1.row(randOrder.at(p)) + hessianXdirection.row(randOrder.at(p));
      d2p_k = Hessian.row(randOrder.at(p)).col(randOrder.at(p));
      
      // if the parameter is regularized:
      if(weights.at(randOrder.at(p)) != 0){
        newParameter = stepDirection.at(randOrder.at(p));
        
        // adjust stepDirection for regularized parameters
        if((dp_k(0,0)-alpha*lambda*weights.row(randOrder.at(p))(0,0)) >= (d2p_k(0,0)*(newParameter(0,0)+stepDirection.row(randOrder.at(p))(0,0)))){
          // condition 1
          z_j = -(dp_k-alpha*lambda*weights.row(randOrder.at(p)))/(d2p_k);
          z.row(randOrder.at(p)) = z_j(0,0);
          stepDirection.row(randOrder.at(p)) += z_j(0,0);
        }else if((dp_k(0,0)+alpha*lambda*weights.row(randOrder.at(p))(0,0))<=(d2p_k(0,0)*(newParameter(0,0)+stepDirection.row(randOrder.at(p))(0,0)))){
          // condition 2
          z_j = -(dp_k+alpha*lambda*weights.row(randOrder.at(p)))/(d2p_k);
          z.row(randOrder.at(p)) = z_j(0,0);
          stepDirection.row(randOrder.at(p)) = stepDirection.row(randOrder.at(p)) + z_j(0,0);
        }else{
          // condition 3
          z_j = -(newParameter+stepDirection.row(randOrder.at(p)));
          z.row(randOrder.at(p)) = z_j(0,0);
          stepDirection.row(randOrder.at(p)) = stepDirection.row(randOrder.at(p))+z_j(0,0);
        }
      }else{
        // if not regularized: coordinate descent with newton direction
        z_j = -dp_k/d2p_k;
        z.row(randOrder.at(p)) = z_j(0,0);
        stepDirection.row(randOrder.at(p)) = stepDirection.row(randOrder.at(p))+z_j(0,0);
      }
    }
    
    // check inner stopping criterion:
    HessTimesZ = HessDiag*pow(z,2);
    
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


inline arma::colvec glmnetLineSearch(
    model& model_,
    penaltyLASSO& penalty_,
    penaltyRidge& smoothPenalty_, 
                             const arma::colvec weights,
                             const double lambda,
                             const double alpha,
                             const arma::colvec parameters_kMinus1, 
                             const double fit_kMinus1, 
                             const arma::colvec gradients_kMinus1,
                             const arma::mat Hessian_kMinus1,
                             const arma::colvec direction,
                             const double stepSize, 
                             const double sigma,
                             const double gamma, 
                             const int maxIterLine){
  arma::colvec newGradients(gradients_kMinus1.n_rows);
  newGradients.fill(arma::datum::nan);
  
  // get penalized M2LL for step size 0:
  penalty_.getValue(parameters_kMinus1);
  pen_0 <- linr::getPenaltyValue(parameters = oldParameters,
                                 lambda = lambda,
                                 alpha = alpha,
                                 regularizedParameterLabels = regularizedParameterLabels,
                                 adaptiveLassoWeights = adaptiveLassoWeights)
  f_0 <- oldM2LL + pen_0
  
  i <- 0
  stepSizeInit <- stepSize
  if(stepSizeInit >= 1) stepSizeInit <- .9
# sometimes the optimizer gets stuck in one specific location;
# here, it can help to try different step sizes
  stepSizeInit <- runif(1,.1,stepSizeInit)
    while(TRUE){
      stepSize <- stepSizeInit^i
      newParameters <- oldParameters+stepSize*direction
      
# compute new fitfunction value
      newM2LL <- try(linr:::fit(linr:::setParameters(SEM = SEM, 
                                                     labels = names(newParameters), 
                                                     values = newParameters, 
                                                     raw = TRUE))$m2LL,
                                                     silent = TRUE)
        if(!is.finite(newM2LL)){
          i <- i+1
          next
        }
        if(is(newM2LL, "try-error")){
          stop("Error in line search: Could not compute fit at current location.")
        }
        
        
# compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
        p_new <- linr::getPenaltyValue(parameters = newParameters,
                                       lambda = lambda,
                                       alpha = alpha,
                                       regularizedParameterLabels = regularizedParameterLabels,
                                       adaptiveLassoWeights = adaptiveLassoWeights)
          f_new <- newM2LL + p_new
          
# test line search criterion
          lineCriterion <- f_new/N - f_0/N <= sig*stepSize*(t(oldGradients/N)%*%direction + 
            gam*t(direction)%*%(oldHessian/N)%*%direction + 
            p_new/N - pen_0/N)
# the likelihoods and the hessian depend on the sample size; re-scaling as otherwise
# the convergence criterion is much more difficult to achieve in larger samples.
# Similarly, getPenaltyValue scales the penalty with the number of subjects and should therefore be re-scaled
# (scaling penalty values with N is handled by adapting lambda and therefore
# not directly visible in the call above) 
            
            if(lineCriterion){
# check if covariance matrix is positive definite
              if(any(eigen(SEM$S, only.values = TRUE)$values < 0)){
                i <- i+1
                next
              }
              if(any(eigen(SEM$impliedCovariance, only.values = TRUE)$values < 0)){
                i <- i+1
                next
              }
# check if gradients can be computed at the new location; this can often cause issues
              newGradients <-  try(linr:::getGradients(SEM, raw = TRUE), silent = TRUE)
                
                if(is(newGradients, "try-error")) {
                  car("Was error\n")
                  i <- i+1
                  next
                }
                break
            }
            i <- i+1
            if(i >= maxIterLine){
              break
            }
    }
    if(is.null(newGradients) || is(newGradients, "try-error")){
      newGradients <- try(linr:::getGradients(SEM, raw = TRUE), silent = TRUE)
    }
    return(
      list("stepSize" = stepSize,
           "newParameters" = newParameters,
           "newGradients" = newGradients)
    )
}


inline fitResults glmnet(model& model_, 
                         arma::rowvec startingValues,
                         penaltyLASSO& penalty_,
                         penaltyRidge& smoothPenalty_, 
                         const tuningParametersEnet tuningParameters, // tuning parameters are of type T
                         const controlGLMNET& control_){
  if(control_.verbose != 0) {
    Rcpp::Rcout << "Optimizing with glmnet.\n" <<
      std::endl;
  }
  
  // we rely heavily on parameter labels because we want to make sure that we 
  // are changing the correct values
  Rcpp::StringVector parameterLabels = startingValues.names();
  
  // prepare parameter vectors
  arma::rowvec parameters_k = Rcpp::clone(startingValues), 
    parameters_kMinus1 = Rcpp::clone(startingValues);
  // the following elements will be required to judge the breaking condition
  arma::rowvec parameterChange(startingValues.length());
  arma::rowvec armaGradients(startingValues.length());
  arma::colvec direction;
  arma::rowvec gradientChange(startingValues.length()); // necessary for Barzilai Borwein
  arma::mat quadr, parchTimeGrad;
  
  // prepare fit elements
  double fit_k = model_.fit(startingValues), 
    fit_kMinus1 = model_.fit(startingValues),
    penalty_k = 0.0;
  double penalizedFit_k, penalizedFit_kMinus1;
  arma::rowvec gradients_k, gradients_kMinus1;
  
  double ridgePenalty = 0.0;
  
  penalizedFit_k = fit_k + 
    penalty_.getValue(parameters_k, tuningParameters) + // lasso penalty part
    smoothPenalty_.getValue(parameters_k, tuningParameters); // ridge penalty part
  penalizedFit_kMinus1 = fit_kMinus1 + 
    penalty_.getValue(parameters_kMinus1, tuningParameters) + // lasso penalty part
    smoothPenalty_.getValue(parameters_kMinus1, tuningParameters); // ridge penalty part
  
  // the following vector will save the fits of all iterations:
  arma::rowvec fits(control_.maxIterOut+1);
  fits.fill(NA_REAL);
  fits.at(0) = penalizedFit_kMinus1;
  
  // prepare gradient elements 
  // NOTE: We combine the gradients of the smooth functions (the log-Likelihood)
  // of the model and the smooth penalty function (e.g., ridge)
  gradients_k = model_.gradients(parameters_k) +
    smoothPenalty_.getGradients(parameters_k, tuningParameters); // ridge part
  gradients_kMinus1 = model_.gradients(parameters_kMinus1) +
    smoothPenalty_.getGradients(parameters_kMinus1, tuningParameters); // ridge part
  
  // prepare Hessian elements
  arma::mat Hessian_k(startingValues.length(), startingValues.length()), 
  Hessian_kMinus1(startingValues.length(), startingValues.length());
  
  Hessian_kMinus1 = control_.initialHessian;
  
  // breaking flags
  bool breakInner = false, // if true, the inner iteration is exited
    breakOuter = false; // if true, the outer iteration is exited
  
  // initialize step size
  double stepSize = control_.stepSize;
  
  // outer iteration
  for(int outer_iteration = 0; outer_iteration < control_.maxIterOut; outer_iteration ++){    
    if(control_.verbose == -99) Rcpp::Rcout << "Outer iteration " << outer_iteration + 1 << std::endl;
    
    // the gradients will be used by the inner iteration to compute the new 
    // parameters
    gradients_kMinus1 = model_.gradients(parameters_kMinus1) +
      smoothPenalty_.getGradients(parameters_kMinus1, tuningParameters); // ridge part
    armaGradients = Rcpp::as<arma::rowvec>(gradients_kMinus1); // needed in convergence criterion
    
    if(control_.verbose == -99) Rcpp::Rcout << "gradients_kMinus1 : " << gradients_kMinus1 << std::endl;
    
    direction = glmnetInner(parameters_kMinus1,
                            Rcpp::as<arma::colvec>(gradients_kMinus1),
                            Hessian_kMinus1,
                            tuningParameters.lambda,
                            tuningParameters.alpha,
                            Rcpp::as<arma::colvec>(tuningParameters.weights),
                            control_.maxIterIn,
                            control_.breakInner);
    
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
        parchTimeGrad = parameterChange*arma::trans(armaGradients); // can be 
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
      
      if(control_.verbose == -99){
        Rcpp::Rcout << "gradients_k\n" << gradients_k << std::endl;
        gradients_k = model_.gradients(parameters_k) + 
          smoothPenalty_.getGradients(parameters_k, tuningParameters);
        Rcpp::Rcout << "Should be identical to\n" << gradients_k << std::endl;
      }
      
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
    penalizedFit_kMinus1 = penalizedFit_k;
    parameters_kMinus1 = Rcpp::clone(parameters_k);
    gradients_kMinus1 = Rcpp::clone(gradients_k);
    
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