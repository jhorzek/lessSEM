#include <RcppArmadillo.h>
#include "implied.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeImpliedCovarianceFull
//' 
//' computes the implied covariance matrix including latent vaiables of a SEM using RAM notation
//' @param Amatrix matrix with directed effects
//' @param Smatrix matrix with undirected effects
//' @param IminusAInverse (I-Amatrix)^(-1)
//' @returns matrix with implied covariances
arma::mat computeImpliedCovarianceFull(const arma::mat& Amatrix, 
                                       const arma::mat& Smatrix,
                                       const arma::mat& IminusAInverse){
  
  return(IminusAInverse*Smatrix*arma::trans(IminusAInverse));
  
}

//' computeImpliedCovariance
//' 
//' computes the implied covariance matrix of a SEM using RAM notation
//' @param Fmatrix filter matrix
//' @param impliedCovarianceFull implied covariance matrix including latent variables
//' @returns matrix with implied covariances
arma::mat computeImpliedCovariance(const arma::mat& Fmatrix, 
                                   const arma::mat& impliedCovarianceFull){
  
  arma::mat impliedCovariance = Fmatrix * impliedCovarianceFull * arma::trans(Fmatrix);
  
  return(impliedCovariance);
}


//' computeImpliedMeansFull
//' 
//' computes the implied means vector of a SEM including the latent variables using RAM notation
//' @param Amatrix matrix with directed effects
//' @param Mvector vector with means
//' @param IminusAInverse (I-Amatrix)^(-1)
//' @returns matrix with implied means
arma::colvec computeImpliedMeansFull(const arma::mat& Amatrix, 
                                  const arma::colvec& Mvector,
                                  const arma::mat& IminusAInverse){
  
  arma::colvec impliedMeans = IminusAInverse*Mvector;
  
  return(impliedMeans);
}

//' computeImpliedMeans
//' 
//' computes the implied means vector of a SEM using RAM notation
//' @param Fmatrix filter matrix
//' @param impliedMeansFull implied means vector including latent variables
//' @returns matrix with implied means
arma::colvec computeImpliedMeans(const arma::mat& Fmatrix, 
                              const arma::mat& impliedMeansFull){
  
  arma::colvec impliedMeans = Fmatrix*impliedMeansFull;
  
  return(impliedMeans);
}

