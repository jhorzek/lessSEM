#include "multivariateNormal.h"

double m2LLMultiVariateNormal(const arma::colvec& observed,
                              const arma::colvec& impliedMeans, 
                              const arma::mat& impliedCovariance,
                              const arma::mat& impliedCovarianceInverse,
                              const double logDetImpliedCovariance){
  
  double m2LL;
  double klog2pi = observed.n_elem*std::log(2*M_PI);
  arma::colvec meanDiff = observed - impliedMeans;
  arma::mat dist = arma::trans(meanDiff)*impliedCovarianceInverse*meanDiff;
  m2LL = klog2pi +
    logDetImpliedCovariance +
    dist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  return(m2LL);
  
}

double m2LLMultiVariateNormalDerivative(
    const std::string& location,
    const arma::colvec& observed,
    const arma::mat& impliedMeans,
    const arma::colvec& impliedMeansDerivative,
    const arma::mat& impliedCovariance, 
    const arma::mat& impliedCovarianceInverse,
    const arma::mat& impliedCovarianceDerivative){
  
  const arma::colvec meanDiff = observed - impliedMeans;
  
  if(location.compare("Mvector") == 0){
    // parameter is in means vector
    
    // \frac{\partial}{\partial \theta_j}(\pmb x - \pmb \mu(\pmb\theta))^T\pmb\Sigma(\pmb\theta)^{-1}(\pmb x - \pmb \mu(\pmb\theta)) :
    const double element_2 = as_scalar(2.0*arma::trans(- impliedMeansDerivative)*
                                       impliedCovarianceInverse * (meanDiff));
    
    return(element_2);
    
  }
  
  
  // parameter is in covariance matrix
  const arma::mat implInvXimplDer = impliedCovarianceInverse * 
    impliedCovarianceDerivative;
  // \frac{\partial}{\partial \theta_j}\ln(|\pmb\Sigma(\pmb\theta)|) : 
  const double element_1 = arma::trace(implInvXimplDer);
  
  // \frac{\partial}{\partial \theta_j}(\pmb x - \pmb \mu(\pmb\theta))^T\pmb\Sigma(\pmb\theta)^{-1}(\pmb x - \pmb \mu(\pmb\theta)) :
  
  if(location.compare("Smatrix") == 0){
    const double element_2 = arma::as_scalar(
      arma::trans(meanDiff)*(-implInvXimplDer) * impliedCovarianceInverse * meanDiff);
    
    return(element_1 + element_2);
  }
  
  if(location.compare("Amatrix") == 0){
    
    const double element_2 = arma::as_scalar(arma::trans(meanDiff)*(
      2.0 * impliedCovarianceInverse * (-impliedMeansDerivative) +
      (-implInvXimplDer) * impliedCovarianceInverse * meanDiff));
    
    return(element_1 + element_2);
    
  }
  
  Rcpp::stop("Unknown parameter location.")
}

}

double m2LLGroupMultiVariateNormal(double N,
                                   const arma::colvec& observedMeans, 
                                   const arma::colvec& impliedMeans, 
                                   const arma::mat& observedCov,
                                   const arma::mat& impliedCovariance,
                                   const arma::mat& impliedCovarianceInverse,
                                   const double logDetImpliedCovariance){
  
  const double Nklog2pi = N*observedMeans.n_elem*std::log(2*M_PI);
  const double NtrSSigma = N*arma::trace(impliedCovarianceInverse * observedCov.t());
  const arma::colvec meanDiff = observedMeans - impliedMeans;
  const arma::mat Ndist = N*arma::trans(meanDiff)*impliedCovarianceInverse * meanDiff;
  
  return(Nklog2pi +
         N*logDetImpliedCovariance +
         NtrSSigma+
         Ndist(0,0)
  ); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  
}

double m2LLGroupMultiVariateNormalDerivative(
    const std::string& location,
    double N,
    const arma::colvec& observedMeans, 
    const arma::colvec& impliedMeans, 
    const arma::colvec& impliedMeansDerivative,
    const arma::mat& observedCov,
    const arma::mat& impliedCovariance,
    const arma::mat& impliedCovarianceInverse,
    const arma::mat& impliedCovarianceDerivative
){
  
  const arma::colvec meanDiff = observedMeans - impliedMeans;
  
  if(location.compare("Mvector") == 0){
    // parameter is in means vector
    
    // \frac{\partial}{\partial \theta_j}(\pmb x - \pmb \mu(\pmb\theta))^T\pmb\Sigma(\pmb\theta)^{-1}(\pmb x - \pmb \mu(\pmb\theta)) :
    double element_3 = as_scalar(N*2.0*arma::trans(-impliedMeansDerivative)*
                                 impliedCovarianceInverse * (meanDiff));
    return(element_3);
    
  }
  
  // parameter is in covariance matrix
  const arma::mat implInvXimplDer = impliedCovarianceInverse * impliedCovarianceDerivative;
  // \frac{\partial}{\partial \theta_j}\ln(|\pmb\Sigma(\pmb\theta)|) : 
  double element_1 = N*arma::trace(implInvXimplDer);
  const double element_2 = N*arma::trace(-observedCov * implInvXimplDer * impliedCovarianceInverse);
  
  if(location.compare("Smatrix") == 0){
    
    // \frac{\partial}{\partial \theta_j}(\pmb x - \pmb \mu(\pmb\theta))^T\pmb\Sigma(\pmb\theta)^{-1}(\pmb x - \pmb \mu(\pmb\theta)) :
    const double element_3 = arma::as_scalar(
      N*
        arma::trans(meanDiff) * (-implInvXimplDer) * impliedCovarianceInverse * meanDiff);
    
    return(element_1 + element_2 + element_3);
  }
  
  if(location.compare("Amatrix") == 0){
    
    // \frac{\partial}{\partial \theta_j}(\pmb x - \pmb \mu(\pmb\theta))^T\pmb\Sigma(\pmb\theta)^{-1}(\pmb x - \pmb \mu(\pmb\theta)) :
    const double element_3 = arma::as_scalar(N*(
      2.0 * arma::trans(meanDiff) * impliedCovarianceInverse * (-impliedMeansDerivative) +
      arma::trans(meanDiff) * (-implInvXimplDer) * impliedCovarianceInverse * meanDiff));
    
    return(element_1 + element_2 + element_3);
    
  }
  
  Rcpp::stop("Unknown parameter location.")
    
}