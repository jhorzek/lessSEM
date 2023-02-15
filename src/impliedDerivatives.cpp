#include "impliedDerivatives.h"

arma::mat impliedCovarianceDerivative(const std::string& location, 
                                      const arma::mat& impliedCovariance, 
                                      const arma::mat& impliedCovarianceFull,
                                      const arma::mat& Fmatrix,
                                      const arma::mat& IminusAInverse,
                                      const arma::mat& derivativeElement){
  
  if(location.compare("Amatrix") == 0){
    arma::mat temp = Fmatrix*IminusAInverse*derivativeElement;
    
    // there is probably some more speed to be found here:
    return(
      temp * impliedCovarianceFull * arma::trans(Fmatrix) +
        Fmatrix * impliedCovarianceFull * arma::trans(temp)
    );
  }
  
  if(location.compare("Smatrix") == 0){
    return(Fmatrix * IminusAInverse * derivativeElement * arma::trans(Fmatrix*IminusAInverse)); 
  }
  
  if(location.compare("Mvector") == 0){
    return(arma::mat(impliedCovariance.n_rows, impliedCovariance.n_cols, arma::fill::zeros));
  }
  
  Rcpp::Rcout << location << std::endl;
  Rcpp::stop("Unknown parameter location");
}

arma::mat impliedMeansDerivative(const std::string& location, 
                                 const arma::mat& impliedMeans, 
                                 const arma::mat& impliedMeansFull,
                                 const arma::mat& Fmatrix,
                                 const arma::mat& IminusAInverse,
                                 const arma::mat& derivativeElement){
  
  if(location.compare("Amatrix") == 0){
    return(-1.0 * Fmatrix * IminusAInverse * ( arma::eye(derivativeElement.n_rows, derivativeElement.n_cols) - derivativeElement ) * impliedMeansFull);
  }
  
  if(location.compare("Smatrix") == 0){
    return(arma::mat(impliedMeans.n_rows, impliedMeans.n_cols, arma::fill::zeros));
  }
  
  if(location.compare("Mvector") == 0){
    return( Fmatrix * IminusAInverse * derivativeElement );
  }
  
  Rcpp::stop("Unknown parameter location");
}