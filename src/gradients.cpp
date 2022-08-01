#include <RcppArmadillo.h>
#include "SEM.h"
#include "scores.h"
#include "gradients.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::rowvec gradientsByGroup(const SEMCpp& SEM, bool raw){

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
  
  // define some elements which will be used repeatedly:
  arma::mat IminusAInverse = arma::inv(arma::eye(arma::size(Amatrix)) - Amatrix);

  std::vector<arma::mat> impliedCovInverses(numberOfMissingnessPatterns);
  for(int m = 0; m < numberOfMissingnessPatterns; m++){
      impliedCovInverses.at(m) = arma::inv_sympd(impliedCovariance(dataSubsets.at(m).notMissing, dataSubsets.at(m).notMissing));
  }
  
  // For each parameter, compute the derivative of the expected covariance
  std::vector<arma::mat>  derivativesOfCovariance = computeimpliedCovarianceDerivatives(SEM,
                                                                                        IminusAInverse,
                                                                                        raw);
  
  // log-det-Sigma-derivative is independent of the person,
  // but depends on the missing structure!
  arma::mat logDetSigmas = computeLogDetSigmas(SEM,
                                               IminusAInverse,
                                               derivativesOfCovariance,
                                               impliedCovInverses,
                                               raw);
  // compute group gradients
  
  
  // individual in row, parameter in column
  arma::mat groupGradients(SEM.data.nGroups, uniqueParameterLabels.size(), arma::fill::zeros);
  arma::rowvec currentGradients(uniqueParameterLabels.size(), arma::fill::zeros);
  arma::colvec dataWithoutMissings; // used in N=1 case
  
  for(int gr = 0; gr < SEM.data.nGroups; gr++){
    
    if(SEM.data.dataSubsets.at(gr).N == 1){
      groupGradients.row(gr) = computeSingleSubjectGradient_Internal(SEM.data.dataSubsets.at(gr).rawData, 
                         SEM,
                         IminusAInverse,
                         dataSubsets.at(gr).notMissing,
                         impliedCovInverses.at(gr),
                         arma::trans(logDetSigmas.row(gr)), 
                         derivativesOfCovariance);
    }
    
    if(SEM.data.dataSubsets.at(gr).N > 1){
      groupGradients.row(gr) = arma::trans(computeSubgroupGradient_Internal(
        SEM.data.dataSubsets.at(gr).N,
        SEM.data.dataSubsets.at(gr).means,
        SEM.data.dataSubsets.at(gr).covariance,
        SEM,
        IminusAInverse,
        dataSubsets.at(gr).notMissing,
        gr,
        impliedCovInverses.at(gr),
        arma::trans(logDetSigmas.row(gr)), 
        derivativesOfCovariance)
      );
    }
    
    currentGradients += groupGradients.row(gr);
    
  }
  
  return(currentGradients);
  
}


arma::mat computeSingleSubjectGradient_Internal(const arma::mat& data_i, // The function also allows
                                                // for a short cut if multiple individuals have the same
                                                // missingness pattern; Then persons are assumed to be in
                                                // the rows, variables in columns.
                                                const SEMCpp& SEM,
                                                const arma::mat& IminusAInverse,
                                                const arma::uvec isObserved,
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
  
  arma::mat individualGradients(data_i.n_rows,uniqueParameterLabels.size(), arma::fill::zeros);
  arma::rowvec tempVec;
  arma::mat tempMat;
  arma::mat dataNoNAMinusImplied = data_i.cols(isObserved);
  dataNoNAMinusImplied.each_row() -= arma::trans(impliedMeans(isObserved));
  
  for(int p = 0; p < uniqueParameterLabels.size(); p++){
    
    if(parameterLocations.at(p).compare("Mvector") == 0){
      
      tempVec = arma::trans(-Fmatrix*IminusAInverse*positionInLocation.at(p));
      individualGradients.col(p) = arma::trans(2.0*arma::trans(tempVec(isObserved))*impliedCovInverse*
        arma::trans(dataNoNAMinusImplied));
      
      continue;
    }
    if(parameterLocations.at(p).compare("Smatrix") == 0){
      
      tempMat = logDetSigmas(p) + dataNoNAMinusImplied*
        (-impliedCovInverse)*
        derivativesOfCovariance.at(p)(isObserved, isObserved)*
        impliedCovInverse*
        arma::trans(dataNoNAMinusImplied);
      individualGradients.col(p) = tempMat.diag();
      
      continue;
    }
    if(parameterLocations.at(p).compare("Amatrix") == 0){
      
      tempVec = arma::trans(-Fmatrix*IminusAInverse*positionInLocation.at(p)*IminusAInverse*Mvector);
      tempMat = dataNoNAMinusImplied*
        (-impliedCovInverse)*
        derivativesOfCovariance.at(p)(isObserved, isObserved)*
        impliedCovInverse*
        arma::trans(dataNoNAMinusImplied);
      
      individualGradients.col(p) =  logDetSigmas(p) + arma::trans(2.0*arma::trans(tempVec(isObserved))*
        impliedCovInverse*arma::trans(dataNoNAMinusImplied)) + tempMat.diag();
      
      continue;
    }
    Rcpp::stop("Encountered parameter in unknown location");
  }
  
  return(individualGradients);
}


arma::colvec computeSubgroupGradient_Internal(
    const int N,
    const arma::colvec& means, 
    const arma::mat& covariance,
    const SEMCpp& SEM,
    const arma::mat& IminusAInverse,
    const arma::uvec isObserved,
    const int groupInSubset,
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
  
  arma::colvec groupGradient(uniqueParameterLabels.size(), arma::fill::zeros);
  arma::rowvec tempVec;
  arma::mat tempMat, tempMat2;
  arma::colvec meanDiff = (means - impliedMeans(isObserved));
  
  for(int p = 0; p < uniqueParameterLabels.size(); p++){
    
    if(parameterLocations.at(p).compare("Mvector") == 0){
      
      tempVec = arma::trans(-Fmatrix*IminusAInverse*positionInLocation.at(p));
      tempMat = arma::trans(tempVec(isObserved))*
        impliedCovInverse*
        (meanDiff);
      groupGradient(p) = N*2.0*tempMat(0,0);
      continue;
    }
    if(parameterLocations.at(p).compare("Smatrix") == 0){
      
      tempMat2 = (-impliedCovInverse)*
        derivativesOfCovariance.at(p)(isObserved, isObserved)*
        impliedCovInverse;
      tempMat = N*logDetSigmas(p) +
        N*arma::trace(covariance*tempMat2)+
        N*arma::trans(meanDiff)*
        tempMat2*
        (meanDiff);
      groupGradient(p) = tempMat(0,0);
      continue;
    }
    if(parameterLocations.at(p).compare("Amatrix") == 0){
      
      tempMat2 = (-impliedCovInverse)*
        derivativesOfCovariance.at(p)(isObserved, isObserved)*
        impliedCovInverse;
      
      tempVec = arma::trans(-Fmatrix*IminusAInverse*positionInLocation.at(p)*IminusAInverse*Mvector);
      tempVec = arma::trans(tempVec(isObserved));
      
      tempMat = N*logDetSigmas(p) +
        N*2*tempVec*impliedCovInverse*
        (meanDiff) +
        N*arma::trace(covariance*tempMat2)+
        N*arma::trans(meanDiff)*
        tempMat2*
        (meanDiff);
      
      groupGradient(p) = tempMat(0,0);
      continue;
    }
    Rcpp::stop("Encountered parameter in unknown location");
  }
  return(groupGradient);
}
