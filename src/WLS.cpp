#include "WLS.h"

double WLS(const arma::mat& weightsInverse,
           const arma::mat& observedCov,
           const arma::colvec& observedMeans,
           const arma::colvec& impliedMeans, 
           const arma::mat& impliedCovariance){
  
  arma::colvec diffMeans = observedMeans - impliedMeans;
  arma::colvec diffCov = arma::vectorise(observedCov - impliedCovariance);
  
  arma::colvec diff = arma::join_rows(diffMeans, diffCov);
  
  return(
    arma::as_scalar(diff.t() * weightsInverse * diff)
  );
  
}

double WLSDerivative(const arma::mat& weightsInverse,
                     const arma::mat& observedCov,
                     const arma::colvec& observedMeans,
                     const arma::colvec& impliedMeans, 
                     const arma::mat& impliedCovariance,
                     const arma::colvec& impliedMeansDerivative,
                     const arma::mat& impliedCovarianceDerivative){
  
  arma::colvec diffMeans = observedMeans - impliedMeans;
  arma::colvec diffCov = arma::vectorise(observedCov - impliedCovariance);
  
  arma::colvec diff = arma::join_rows(diffCov, diffMeans);
  
  arma::colvec impl = arma::join_rows(impliedMeansDerivative, 
                                      arma::vectorise(impliedCovarianceDerivative));
  
  return(arma::as_scalar(
      arma::trans((weightsInverse + arma::trans(weightsInverse)) * diff) * impl
  ));
}
