#include <RcppArmadillo.h>
#include "config.hpp"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the -2log likelihood in a RAM model for a  single person given the number of non-missing variables (nObservedVariables), x (rawData), filtered impliedMeans, and filtered impliedCovariance.

// [[Rcpp::export]]
double computeIndividualM2LL(const int nObservedVariables, const arma::colvec& rawData, const arma::colvec& impliedMeans, const arma::mat& impliedCovariance){
  double m2LL;
  double klog2pi = nObservedVariables*std::log(2*M_PI);
  double logDetExpCov = std::log(arma::det(impliedCovariance));
  arma::mat dist = arma::trans(rawData - impliedMeans)*arma::inv(impliedCovariance)*(rawData - impliedMeans);
  m2LL = klog2pi +
    logDetExpCov +
    dist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  return(m2LL);
}

// Computes the -2log likelihood in a RAM model for a  group of people with identical missing structure given the sample size, number of non-missing variables (nObservedVariables),
// observedMeans within this sample, observedCov within this sample, filtered impliedMeans, and filtered impliedCovariance.
// The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05

// [[Rcpp::export]]
double computeGroupM2LL(const int sampleSize, const int nObservedVariables, const arma::colvec& observedMeans, const arma::mat& observedCov,
                        const arma::colvec& impliedMeans, const arma::mat& impliedCovariance){
  double m2LL;
  double Nklog2pi = sampleSize*nObservedVariables*std::log(2*M_PI);
  double NlogDetExpCov = sampleSize*std::log(arma::det(impliedCovariance));
  double NtrSSigma = sampleSize*arma::trace(observedCov * arma::inv(impliedCovariance));
  arma::mat Ndist = sampleSize*arma::trans(observedMeans - impliedMeans)*arma::inv(impliedCovariance)*(observedMeans - impliedMeans);

  m2LL = Nklog2pi +
    NlogDetExpCov +
    NtrSSigma+
    Ndist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible

  return(m2LL);
}
