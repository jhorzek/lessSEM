#include <RcppArmadillo.h>
#include "fit.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeIndividualM2LL
//' 
//' Computes the -2log likelihood in a SEM for a single person 
//' @param nObservedVariables number of non-missing variables 
//' @param rawData observed data for this person
//' @param impliedMeans implied means of the SEM
//' @param impliedCovariance implied covariance of the SEM
//' @returns -2 log-Lilekelihood
//' 
// [[Rcpp::export]]
double computeIndividualM2LL(const int nObservedVariables, 
                             const arma::colvec& rawData, 
                             const arma::colvec& impliedMeans, 
                             const arma::mat& impliedCovariance){
  double m2LL;
  double klog2pi = nObservedVariables*std::log(2*M_PI);
  double logDetExpCov = std::log(arma::det(impliedCovariance));
  arma::mat dist = arma::trans(rawData - impliedMeans)*arma::inv(impliedCovariance)*(rawData - impliedMeans);
  m2LL = klog2pi +
    logDetExpCov +
    dist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  return(m2LL);
}


//' computeGroupM2LL
//' Computes the -2log likelihood in a RAM model for a  group of people with identical 
//' missing structure. The idea is based on OpenMx (Boker, S. M., Neale, M. C., Maes, H. H., 
//' Wilde, M. J., Spiegel, M., Brick, T. R., Estabrook, R., Bates, T. C., & Mehta, P. (2022). 
//' OpenMx: Extended structural equation modelling (2.20.0). https://cran.r-project.org/web/packages/OpenMx/OpenMx.pdf)
//' @param sampleSize number of persons
//' @param nObservedVariables number of non-missing variables
//' @param observedMeans vector with means of the observed variables
//' @param observedCov matrix with (co-)variances of the observed variables
//' @param impliedMeans implied means of the SEM
//' @param impliedCovariance implied covariance of the SEM
//' @returns -2 log-Lilekelihood
//' 
// [[Rcpp::export]]
double computeGroupM2LL(const int sampleSize, 
                        const int nObservedVariables,
                        const arma::colvec& observedMeans, 
                        const arma::mat& observedCov,
                        const arma::colvec& impliedMeans,
                        const arma::mat& impliedCovariance){
  double m2LL;
  arma::mat inv = arma::inv(impliedCovariance);
  double Nklog2pi = sampleSize*nObservedVariables*std::log(2*M_PI);
  double NlogDetExpCov = sampleSize*std::log(arma::det(impliedCovariance));
  double NtrSSigma = sampleSize*arma::trace(observedCov * inv);
  arma::mat Ndist = sampleSize*arma::trans(observedMeans - impliedMeans)*inv*(observedMeans - impliedMeans);

  m2LL = Nklog2pi +
    NlogDetExpCov +
    NtrSSigma+
    Ndist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible

  return(m2LL);
}
