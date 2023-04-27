#include "WLS.h"

// Weighted Least Squares when no means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @return double indicating fit 
double WLS(const arma::mat& weightsInverse,
           const arma::mat& observedCov,
           const arma::mat& impliedCovariance){
  
  arma::colvec diff((observedCov.n_rows*(observedCov.n_cols+1))/2);
  
  unsigned int to = 0;
  for(unsigned int i = 0; i < observedCov.n_rows; i++){
    for(unsigned int j = i; j < observedCov.n_rows; j++){
      
      diff(to) = observedCov(i,j) - impliedCovariance(i,j);
      to++;
    }
  }
  
  return(
    arma::as_scalar(diff.t() * weightsInverse * diff)
  );
  
}

// Weighted Least Squares when means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @parma observedMeans observed means
// @param impliedMeans model implied means
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @return double indicating fit 
double WLS(const arma::mat& weightsInverse,
           const arma::colvec& observedMeans,
           const arma::colvec& impliedMeans, 
           const arma::mat& impliedCovariance,
           const arma::mat& observedCov){
  
  arma::colvec diff((observedCov.n_rows*(observedCov.n_cols+1))/2 + observedCov.n_rows);
  
  unsigned int to = 0;
  for(unsigned int i = 0; i < observedMeans.n_elem; i++){
    diff(to) = observedMeans(i) - impliedMeans(i);
    to++;
  }
  
  for(unsigned int i = 0; i < observedCov.n_rows; i++){
    for(unsigned int j = i; j < observedCov.n_rows; j++){
      diff(to) = observedCov(i,j) - impliedCovariance(i,j);
      to++;
    }
  }
  
  return(
    arma::as_scalar(diff.t() * weightsInverse * diff)
  );
  
}

// WLSDerivative Least Squares when no means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @param impliedCovarianceDerivative derivative of covariance matrix wrt parameter
// @return double indicating derivative of functions wrt individual parameter 
double WLSDerivative(const arma::mat& weightsInverse,
                     const arma::mat& observedCov,
                     const arma::mat& impliedCovariance,
                     const arma::mat& impliedCovarianceDerivative){
  
  arma::rowvec diff((observedCov.n_rows*(observedCov.n_cols+1))/2);
  arma::colvec impl((observedCov.n_rows*(observedCov.n_cols+1))/2);
  
  unsigned int to = 0;
  for(unsigned int i = 0; i < observedCov.n_rows; i++){
    for(unsigned int j = i; j < observedCov.n_rows; j++){
      
      diff(to) = impliedCovariance(i,j) - observedCov(i,j);
      impl(to) = impliedCovarianceDerivative(i,j);
      to++;
    }
  }
  
  return(arma::as_scalar(
      (2.0*diff*weightsInverse) * impl
  ));
}

// WLSDerivative Least Squares when means are estimated
//
// @param weightsInverse matrix with weights from lavaan. First means then covariances
// @parma observedMeans observed means
// @param impliedMeans model implied means
// @param observedCov observed covariance matrix
// @param impliedCovariance implied covariance matrix
// @param impliedCovarianceDerivative derivative of covariance matrix wrt parameter
// @return double indicating derivative of functions wrt individual parameter 
double WLSDerivative(const arma::mat& weightsInverse,
                     const arma::mat& observedCov,
                     const arma::colvec& observedMeans,
                     const arma::colvec& impliedMeans, 
                     const arma::mat& impliedCovariance,
                     const arma::colvec& impliedMeansDerivative,
                     const arma::mat& impliedCovarianceDerivative){
  
  
  arma::rowvec diff((observedCov.n_rows*(observedCov.n_cols+1))/2 + observedCov.n_rows);
  arma::colvec impl((observedCov.n_rows*(observedCov.n_cols+1))/2 + observedCov.n_rows);
  
  unsigned int to = 0;
  
  for(unsigned int i = 0; i < observedCov.n_rows; i++){
    
    diff(to) = impliedMeans(i) - observedMeans(i);
    impl(to) = impliedMeansDerivative(i);
    to++;
  }
  
  for(unsigned int i = 0; i < observedCov.n_rows; i++){
    for(unsigned int j = i; j < observedCov.n_rows; j++){
      
      diff(to) = impliedCovariance(i,j) - observedCov(i,j);
      impl(to) = impliedCovarianceDerivative(i,j);
      to++;
    }
  }
  
  return(arma::as_scalar(
      (2.0*diff*weightsInverse) * impl
  ));
}