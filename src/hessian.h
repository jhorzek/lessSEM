// THE FOLLOWING CODE IS ADAPTED FROM LAVAAN. 
// SEE lavaan:::lav_model_hessian FOR THE IMPLEMENTATION
// BY Yves Rosseel. The code is under GPL (>= 2)

#ifndef HESSIAN_H
#define HESSIAN_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

template<typename T> // T is the type of the model

arma::mat approximateHessian(T& SEM, 
                             Rcpp::StringVector label_,
                             arma::vec value_,
                             bool raw,
                             double eps){
  int nPar = label_.length();
  arma::mat hessian(nPar, nPar, arma::fill::zeros);
  arma::vec stepLeft = value_, 
    twoStepLeft = value_, 
    stepRight = value_, 
    twoStepRight = value_;
  
  arma::rowvec gradientsStepLeft(nPar);
  arma::rowvec gradientsTwoStepLeft(nPar);
  arma::rowvec gradientsStepRight(nPar);
  arma::rowvec gradientsTwoStepRight(nPar);
  
  for(int p = 0; p < nPar; p++) {
    
    stepLeft.at(p) -= eps;
    twoStepLeft.at(p) -= 2*eps;
    stepRight.at(p) += eps;
    twoStepRight.at(p) += 2*eps;
    
    // step left
    SEM.setParameters(label_, stepLeft, raw);
    SEM.fit();
    gradientsStepLeft = SEM.getGradients(raw);
    
    // two step left
    SEM.setParameters(label_, twoStepLeft, raw);
    SEM.fit();
    gradientsTwoStepLeft = SEM.getGradients(raw);
    
    // step right
    SEM.setParameters(label_, stepRight, raw);
    SEM.fit();
    gradientsStepRight = SEM.getGradients(raw);
    
    // two step right
    SEM.setParameters(label_, twoStepRight, raw);
    SEM.fit();
    gradientsTwoStepRight = SEM.getGradients(raw);
    
    // approximate hessian
    hessian.col(p) = arma::trans((gradientsTwoStepLeft - 
      8.0 * gradientsStepLeft + 
      8.0 * gradientsStepRight - 
      gradientsTwoStepRight)/(12.0 * eps));
    
    // reset
    stepLeft.at(p) += eps;
    twoStepLeft.at(p) += 2*eps;
    stepRight.at(p) -= eps;
    twoStepRight.at(p) -= 2*eps;
  }
  // make symmetric
  hessian = (hessian + arma::trans(hessian))/2.0;
  
  // reset SEM
  SEM.setParameters(label_, value_, raw);
  SEM.fit();
  return(hessian);
}

#endif