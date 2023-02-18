#include "impliedDerivatives.h"

arma::mat impliedCovarianceDerivative(const double parameterValue, // only required for variances due to the internal transformation used by lessSEM.
                                      const std::string& location, 
                                      const bool isVariance, // only required for variances due to the internal transformation used by lessSEM.
                                      const bool raw, // only required for variances due to the internal transformation used by lessSEM.
                                      const arma::mat& impliedCovariance, 
                                      const derivPrecompute& precomputedElements,
                                      const arma::mat& Fmatrix,
                                      const arma::mat& IminusAInverse,
                                      const arma::mat& derivativeElement){
  
  if(location.compare("Amatrix") == 0){
    arma::mat temp = precomputedElements.FIminusAInverse*derivativeElement;
    
    // there is probably some more speed to be found here:
    return(
      temp * precomputedElements.impliedCovarianceFulltF +
        precomputedElements.FimpliedCovarianceFull * arma::trans(temp)
    );
  }
  
  if(location.compare("Smatrix") == 0){
    if(isVariance && raw){
      // note: lessSEM uses a log-transform internally. This must be accounted for here.
      // find value of variance
      
      return(precomputedElements.FIminusAInverse * (parameterValue*derivativeElement)*precomputedElements.tFIminusAInverse);
      
    }
    
    return(precomputedElements.FIminusAInverse * derivativeElement * precomputedElements.tFIminusAInverse); 
  }
  
  if(location.compare("Mvector") == 0){
    return(arma::mat(impliedCovariance.n_rows, impliedCovariance.n_cols, arma::fill::zeros));
  }
  
  Rcpp::Rcout << location << std::endl;
  Rcpp::stop("Unknown parameter location");
}

arma::colvec impliedMeansDerivative(const std::string& location, 
                                 const arma::mat& impliedMeans, 
                                 const arma::mat& impliedMeansFull,
                                 const arma::mat& Fmatrix,
                                 const derivPrecompute& precomputedElements,
                                 const arma::mat& IminusAInverse,
                                 const arma::mat& derivativeElement){
  
  if(location.compare("Amatrix") == 0){
    return(-1.0 * precomputedElements.FIminusAInverse * derivativeElement * impliedMeansFull);
  }
  
  if(location.compare("Smatrix") == 0){
    return(arma::mat(impliedMeans.n_rows, impliedMeans.n_cols, arma::fill::zeros));
  }
  
  if(location.compare("Mvector") == 0){
    return( precomputedElements.FIminusAInverse * derivativeElement );
  }
  
  Rcpp::stop("Unknown parameter location");
}