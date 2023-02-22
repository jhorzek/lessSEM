#ifndef GRADIENTS_H
#define GRADIENTS_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "SEM.h"
#include "impliedDerivatives.h"

// [[Rcpp :: depends ( RcppArmadillo , RcppParallel)]]

class initializeGradients: public RcppParallel::Worker{
public:
  SEMCpp& SEM;
  const bool raw;
  int nParameters;
  derivPrecompute precomputedElements;
  
  initializeGradients(SEMCpp& SEM, const bool raw): SEM(SEM), raw(raw){
    
    nParameters = SEM.derivElements.uniqueLabels.size();
    
    if(!SEM.detivativesInitialized){
      // we have to initialize the subgroup gradients first
      
      SEM.impliedCovarianceDerivatives.resize(nParameters);
      SEM.impliedMeansDerivatives.resize(nParameters);
      
      SEM.detivativesInitialized = true;
    }
    
    // We will start by computing the derivatives of the implied covariance
    // matrix and the implied means vector for every parameter. Some elements 
    // are used repeatedly, so it makes sense to just compute them once and
    // re-use them:        
    precomputedElements.FIminusAInverse = SEM.Fmatrix * SEM.IminusAInverse;
    precomputedElements.tFIminusAInverse = arma::trans(precomputedElements.FIminusAInverse);
    precomputedElements.FimpliedCovarianceFull = SEM.Fmatrix * SEM.impliedCovarianceFull;
    precomputedElements.impliedCovarianceFulltF = SEM.impliedCovarianceFull*arma::trans(SEM.Fmatrix);
  }
  
  // compute the gradients in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t p = begin; p < end; p++){
      
      const double parameterValue = SEM.parameterTable.parameterMap.at(
        SEM.derivElements.uniqueLabels.at(p)).value;
      const bool isVariance = SEM.parameterTable.parameterMap.at(
        SEM.derivElements.uniqueLabels.at(p)).isVariance;
      
      SEM.impliedCovarianceDerivatives.at(p) = 
        impliedCovarianceDerivative(
          parameterValue,
          SEM.derivElements.uniqueLocations.at(p),
          isVariance,
          raw,
          SEM.impliedCovariance,
          precomputedElements,
          SEM.Fmatrix,
          SEM.IminusAInverse,
          SEM.derivElements.positionInLocation.at(p)
        );
      
      SEM.impliedMeansDerivatives.at(p) =
        impliedMeansDerivative(SEM.derivElements.uniqueLocations.at(p),
                               SEM.impliedMeans,
                               SEM.impliedMeansFull,
                               SEM.Fmatrix,
                               precomputedElements,
                               SEM.IminusAInverse,
                               SEM.derivElements.positionInLocation.at(p));
    }
  };  
};

class ParallelGradients: public RcppParallel::Worker{
public:
  const SEMCpp& SEM;
  const int mp; // specific missingness pattern -> the group for which the gradients are computed
  const bool raw;
  int nParameters;
  
  arma::rowvec groupGradients;
  
  ParallelGradients(const SEMCpp& SEM, const int mp, const bool raw): SEM(SEM), mp(mp), raw(raw){
    nParameters = SEM.derivElements.uniqueLabels.size();
    groupGradients.resize(nParameters);
    groupGradients.fill(arma::datum::nan);
  }
  
  // compute the gradients in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t p = begin; p < end; p++){
      // compute the gradient for a specific parameter p:
      
      if(SEM.data.dataSubsets.at(mp).N == 1){
        
        groupGradients.col(p) = m2LLMultiVariateNormalDerivative(
          SEM.derivElements.uniqueLocations.at(p),
          SEM.data.dataSubsets.at(mp).rawData(SEM.data.dataSubsets.at(mp).notMissing),
          SEM.subsetImpliedMeans.at(mp),
          SEM.impliedMeansDerivatives.at(p).rows(SEM.data.dataSubsets.at(mp).notMissing),
          SEM.subsetImpliedCovariance.at(mp),
          SEM.subsetImpliedCovarianceInverse.at(mp),
          SEM.impliedCovarianceDerivatives.at(p).submat(SEM.data.dataSubsets.at(mp).notMissing,
                                              SEM.data.dataSubsets.at(mp).notMissing));
      }else{
        
        groupGradients.col(p) = m2LLGroupMultiVariateNormalDerivative(
          SEM.derivElements.uniqueLocations.at(p),
          SEM.data.dataSubsets.at(mp).N,
          SEM.data.dataSubsets.at(mp).means,
          SEM.subsetImpliedMeans.at(mp),
          SEM.impliedMeansDerivatives.at(p).rows(SEM.data.dataSubsets.at(mp).notMissing),
          SEM.data.dataSubsets.at(mp).covariance,
          SEM.subsetImpliedCovariance.at(mp),
          SEM.subsetImpliedCovarianceInverse.at(mp),
          SEM.impliedCovarianceDerivatives.at(p).submat(SEM.data.dataSubsets.at(mp).notMissing,
                                              SEM.data.dataSubsets.at(mp).notMissing));
      }
    }  
  };
};

arma::rowvec gradientsByGroup(const SEMCpp& SEM, bool raw);

#endif