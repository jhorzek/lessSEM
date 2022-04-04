#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "scores.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::colvec computeSingleSubjectGradient_Internal(const arma::colvec& data_i, 
                                                  const SEMCpp& SEM,
                                                  const arma::mat& IminusAInverse,
                                                  const arma::uvec isObserved,
                                                  const int personInMissingGroup,
                                                  const arma::mat& impliedCovInverse,
                                                  const arma::colvec& logDetSigmas, 
                                                  const std::vector<arma::mat>& derivativesOfCovariance){

  // define some convenient references
  const arma::mat& Amatrix = SEM.Amatrix;
  const arma::mat& Smatrix = SEM.Smatrix;
  const arma::mat& Fmatrix = SEM.Fmatrix;
  const arma::colvec& Mvector = SEM.Mvector;
  const arma::mat& impliedCovariance = SEM.impliedCovariance;
  const arma::colvec& impliedMeans = SEM.impliedMeans;
  
  const std::vector<std::string>& uniqueParameterLabels = SEM.derivElements.uniqueLabels;
  const std::vector<std::string>& parameterLocations = SEM.derivElements.uniqueLocations;
  const std::vector<bool>& isVariance = SEM.derivElements.isVariance;
  const std::vector<arma::mat>& positionInLocation = SEM.derivElements.positionInLocation;
  
  const int numberOfMissingnessPatterns = SEM.data.nGroups;
  const std::vector<subset>& dataSubsets = SEM.data.dataSubsets;

  arma::colvec individualGradient(uniqueParameterLabels.size(), arma::fill::zeros);
  arma::rowvec tempVec;
  arma::mat tempMat;

  for(int p = 0; p < uniqueParameterLabels.size(); p++){

    if(parameterLocations.at(p).compare("Mvector") == 0){

      tempVec = arma::trans(-Fmatrix*IminusAInverse*positionInLocation.at(p));
      tempMat = arma::trans(tempVec(isObserved))*
        impliedCovInverse*
        (data_i(isObserved) - impliedMeans(isObserved));
      individualGradient(p) = 2.0*tempMat(0,0);

      continue;
    }
    if(parameterLocations.at(p).compare("Smatrix") == 0){

      tempMat = logDetSigmas(p) + arma::trans(data_i(isObserved) - impliedMeans(isObserved))*
        (-impliedCovInverse)*
        derivativesOfCovariance.at(p)(isObserved, isObserved)*
        impliedCovInverse*
        (data_i(isObserved) - impliedMeans(isObserved));
      individualGradient(p) = tempMat(0,0);
      
      continue;
    }
    if(parameterLocations.at(p).compare("Amatrix") == 0){

      tempVec = arma::trans(-Fmatrix*IminusAInverse*positionInLocation.at(p)*IminusAInverse*Mvector);

      tempMat = logDetSigmas(p) + 2.0*arma::trans(tempVec(isObserved))*
        impliedCovInverse*(data_i(isObserved) - impliedMeans(isObserved))+
        arma::trans(data_i(isObserved) - impliedMeans(isObserved))*
        (-impliedCovInverse)*
        derivativesOfCovariance.at(p)(isObserved, isObserved)*
        impliedCovInverse*
        (data_i(isObserved) - impliedMeans(isObserved));
      individualGradient(p) = tempMat(0,0);
      
      continue;
    }
    Rcpp::stop("Encountered parameter in unknown location");
  }
  return(individualGradient);
}