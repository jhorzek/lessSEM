#include <RcppArmadillo.h>
#include "SEM.h"
#include "scores.h"
#include "multivariateNormal.h"
#include "gradients.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat scores(const SEMCpp& SEM, bool raw){
  
  // define some convenient references
  const std::vector<std::string>& uniqueParameterLabels = SEM.derivElements.uniqueLabels;
  const int numberOfMissingnessPatterns = SEM.data.nGroups;
  const std::vector<subset>& dataSubsets = SEM.data.dataSubsets;
  
  // compute individual gradients
  
  const int N = SEM.rawData.n_rows;

  // individual in row, parameter in column
  arma::mat individualGradients(N, uniqueParameterLabels.size());
  arma::rowvec currentGradients(uniqueParameterLabels.size(), arma::fill::zeros);
  
  arma::uvec rowi(1);
  
  // we can speed up the computation considerably by going over the missingness patterns
  for(unsigned int miss = 0; miss < numberOfMissingnessPatterns; miss++){
    
    // iterate over individuals
    for(arma::uword i = 0; i < dataSubsets.at(miss).N; i++){
      rowi.at(0) = i;
      // iterate over parameters
      for(unsigned int p = 0; p < uniqueParameterLabels.size(); p++){
        
        currentGradients.col(p) = m2LLMultiVariateNormalDerivative(
          SEM.derivElements.uniqueLocations.at(p),
          arma::trans(dataSubsets.at(miss).rawData.submat(rowi, dataSubsets.at(miss).notMissing)),
          SEM.subsetImpliedMeans.at(miss),
          SEM.impliedMeansDerivatives.at(p).rows(dataSubsets.at(miss).notMissing),
          SEM.subsetImpliedCovariance.at(miss),
          SEM.subsetImpliedCovarianceInverse.at(miss),
          SEM.impliedCovarianceDerivatives.at(p).submat(dataSubsets.at(miss).notMissing,
                                              dataSubsets.at(miss).notMissing)
        );
      }

      individualGradients.row(i) = currentGradients;

    }
    
  }

  return(individualGradients);
  
}