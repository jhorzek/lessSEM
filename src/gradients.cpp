#include <RcppArmadillo.h>
#include "SEM.h"
#include "multivariateNormal.h"
#include "scores.h"
#include "gradients.h"
#include "impliedDerivatives.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

void initializeGradients(SEMCpp& SEM, const bool raw){
  
  const int nParameters = SEM.derivElements.uniqueLabels.size();
  
  if(!SEM.detivativesInitialized){
    // we have to initialize the subgroup gradients first
    
    SEM.impliedCovarianceDerivatives.resize(nParameters);
    SEM.impliedMeansDerivatives.resize(nParameters);
    
    SEM.detivativesInitialized = true;
  }
  
  // We will start by computing the derivatives of the implied covariance
  // matrix and the implied means vector for every parameter.
  
  for(int p = 0; p < nParameters; p++){
    
    double parameterValue = SEM.parameterTable.parameterMap.at(
      SEM.derivElements.uniqueLabels.at(p)).value;
    bool isVariance = SEM.parameterTable.parameterMap.at(
      SEM.derivElements.uniqueLabels.at(p)).isVariance;
    
    SEM.impliedCovarianceDerivatives.at(p) = 
      impliedCovarianceDerivative(
        parameterValue,
        SEM.derivElements.uniqueLocations.at(p),
        isVariance,
        raw,
        SEM.impliedCovariance,
        SEM.impliedCovarianceFull,
        SEM.Fmatrix,
        SEM.IminusAInverse,
        SEM.derivElements.positionInLocation.at(p)
      );
    
    SEM.impliedMeansDerivatives.at(p) =
      impliedMeansDerivative(SEM.derivElements.uniqueLocations.at(p),
                             SEM.impliedMeans,
                             SEM.impliedMeansFull,
                             SEM.Fmatrix,
                             SEM.IminusAInverse,
                             SEM.derivElements.positionInLocation.at(p));
  }
  
}

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
          SEM.impliedMeansDerivatives.at(p).rows(dataSubsets.at(mp).notMissing),
          SEM.subsetImpliedCovariance.at(mp),
          SEM.subsetImpliedCovarianceInverse.at(mp),
          SEM.impliedCovarianceDerivatives.at(p).submat(dataSubsets.at(mp).notMissing,
                                              dataSubsets.at(mp).notMissing));
      }else{
        
        gradients.col(p) += m2LLGroupMultiVariateNormalDerivative(
          uniqueLocations.at(p),
          dataSubsets.at(mp).N,
          dataSubsets.at(mp).means,
          SEM.subsetImpliedMeans.at(mp),
          SEM.impliedMeansDerivatives.at(p).rows(dataSubsets.at(mp).notMissing),
          dataSubsets.at(mp).covariance,
          SEM.subsetImpliedCovariance.at(mp),
          SEM.subsetImpliedCovarianceInverse.at(mp),
          SEM.impliedCovarianceDerivatives.at(p).submat(dataSubsets.at(mp).notMissing,
                                              dataSubsets.at(mp).notMissing));
      }
      
    }
    
  }
  
  return(gradients);
  
}