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