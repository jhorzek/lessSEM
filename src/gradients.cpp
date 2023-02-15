#include <RcppArmadillo.h>
#include "SEM.h"
#include "multivariateNormal.h"
#include "scores.h"
#include "gradients.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::rowvec gradientsByGroup(const SEMCpp& SEM, bool raw){
  
  // now that we have initialized all the derivatives of the implied means and
  // covariances, we can start computing the actual gradients
  
  const int numberOfMissingnessPatterns = SEM.data.nGroups;
  const std::vector<std::string>& uniqueParameterLabels = SEM.derivElements.uniqueLabels;
  const std::vector<std::string>& uniqueLocations = SEM.derivElements.uniqueLocations;
  
  const int nParameters = uniqueParameterLabels.size();
  
  arma::rowvec gradients(nParameters, arma::fill::zeros);
  
  const std::vector<subset>& dataSubsets = SEM.data.dataSubsets;
  
  for(int mp = 0; mp < numberOfMissingnessPatterns; mp++){
    
    for(int p = 0; p < nParameters; p++){
      
      if(dataSubsets.at(mp).N == 1){
        gradients.col(p) += m2LLMultiVariateNormalDerivative(
          uniqueLocations.at(p),
          dataSubsets.at(mp).rawData(dataSubsets.at(mp).notMissing),
          SEM.subsetImpliedMeans.at(mp),
          SEM.subsetImpliedMeansDerivatives.at(mp).at(p),
          SEM.subsetImpliedCovariance.at(mp),
          SEM.subsetImpliedCovarianceInverse.at(mp),
          SEM.subsetImpliedCovarianceDerivatives.at(mp).at(p));
      }else{
        gradients.col(p) += m2LLGroupMultiVariateNormalDerivative(
          uniqueLocations.at(p),
          dataSubsets.at(mp).N,
          dataSubsets.at(mp).means,
          SEM.subsetImpliedMeans.at(mp),
          SEM.subsetImpliedMeansDerivatives.at(mp).at(p),
          dataSubsets.at(mp).covariance,
          SEM.subsetImpliedCovariance.at(mp),
          SEM.subsetImpliedCovarianceInverse.at(mp),
          SEM.subsetImpliedCovarianceDerivatives.at(mp).at(p));
      }
      
    }
    
  }
  
  return(gradients);
  
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
