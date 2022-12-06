#ifndef BFGS_H
#define BFGS_H

namespace lessSEM{

// BFGS
//
// computes the BFGS Hessian approximation
//
//
// param parameters_kMinus1 parameters of previous iteration
// param gradients_kMinus1 gradients of previous iteration
// param Hessian_kMinus1 Hessian of previous iteration
// param parameters_k parameters of current iteration
// param gradients_k gradients of current iteration
// param cautious boolean: should the update be skipped if it would result in a non positive definite Hessian?
// param hessianEps controls when the update of the Hessian approximation is skipped
inline arma::mat BFGS(
    const arma::rowvec& parameters_kMinus1, 
    const arma::rowvec& gradients_kMinus1, 
    const arma::mat& Hessian_kMinus1, 
    const arma::rowvec& parameters_k, 
    const arma::rowvec& gradients_k, 
    const bool cautious, 
    const double hessianEps,
    bool verbose){
  
  arma::rowvec y = gradients_k - gradients_kMinus1;
  arma::rowvec d = parameters_k - parameters_kMinus1;
  arma::mat yTimesD = y*arma::trans(d);
  arma::mat Hessian_k = Hessian_kMinus1;
  arma::mat ySquared, dHd;
  bool skipUpdate = false;
  // test if positive definiteness is ensured
  try{
    
    skipUpdate = (yTimesD(0,0) < hessianEps) && cautious;
    
  }catch(...){
    // skip in case of error: return Hessian_kMinus1
    if(verbose) Rcpp::warning("Hessian update skipped.");
    return(Hessian_k);
  }
  
  if(yTimesD(0,0) < 0){
    if(verbose) Rcpp::warning("Hessian update possibly non-positive definite.");
    if(skipUpdate) return(Hessian_k);
  }
  
  // see e.g., Nocedal, J., & Wright, S. J. (2006). Numerical optimization (2nd ed).
  // Springer, p. 537 Equation 18.16
  
  ySquared = arma::trans(y)*y;
  dHd = d*Hessian_kMinus1*arma::trans(d);
  
  Hessian_k = Hessian_kMinus1 - 
    (Hessian_kMinus1*arma::trans(d)*d*Hessian_kMinus1)/dHd(0,0) + 
    ySquared/yTimesD(0,0);
  
  if(!arma::is_finite(Hessian_k)){
    if(verbose) Rcpp::warning("Non-finite Hessian. Returning previous Hessian");
    return(Hessian_kMinus1);
  }
  
  // check for symmetric positive definiteness 
  if(!Hessian_k.is_symmetric()){
    // make symmetric
    double sumElem = arma::accu(arma::pow(Hessian_k - .5*(Hessian_k + arma::trans(Hessian_k)),2));
    if((sumElem > 1) & verbose) Rcpp::warning("Hessian not symmetric");
    Hessian_k = .5*(Hessian_k + arma::trans(Hessian_k));
  }else{
    return(Hessian_k);
  }
  
  // we now know that the matrix is symmetric; lets check again
  // for positive definite
  if(!Hessian_k.is_sympd()){
    // make positive definite
    // see https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
    if(verbose) Rcpp::warning("Hessian not pd");
    arma::vec eigenValues = arma::eig_sym(Hessian_k);
    arma::mat diagMat = arma::eye(Hessian_k.n_rows, 
                                  Hessian_k.n_cols);
    diagMat.fill(0.0);
    diagMat.diag() += -1.1*arma::min(eigenValues);
    Hessian_k = Hessian_k +  diagMat;
    
    // check again...
    if(!Hessian_k.is_sympd()){
      // return non-updated hessian
      if(verbose) Rcpp::warning("Invalid Hessian. Returning previous Hessian");
      return(Hessian_kMinus1);
    }
  }
  
  return(Hessian_k);
  
}

}

#endif