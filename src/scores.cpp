#include <RcppArmadillo.h>
#include "config.hpp"
#include "SEM.hpp"
#include "scores.hpp"
#include "gradients.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat scores(const SEMCpp& SEM, bool raw){
  
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
    arma::mat currentimpliedCovarianceInverse = arma::inv_sympd(impliedCovariance(dataSubsets.at(m).notMissing, dataSubsets.at(m).notMissing));
    
    impliedCovInverses.at(m) = currentimpliedCovarianceInverse;
    
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

  // compute individual gradients
  
  const int N = SEM.rawData.n_rows;

  // individual in row, parameter in column
  arma::mat individualGradients(N, uniqueParameterLabels.size());
  arma::colvec currentGradients;

  for(int i = 0; i < N; i++){
    
    const int personInSubset = SEM.personInSubset.at(i);
    
    currentGradients = computeSingleSubjectGradient_Internal(arma::trans(SEM.rawData.row(i)), 
                                                             SEM,
                                                             IminusAInverse,
                                                             dataSubsets.at(personInSubset).notMissing,
                                                             personInSubset,
                                                             impliedCovInverses.at(personInSubset),
                                                             arma::trans(logDetSigmas.row(personInSubset)), 
                                                             derivativesOfCovariance);
    
    individualGradients.row(i) = arma::trans(currentGradients);
    
  }
  return(individualGradients);
  
}




std::vector<arma::mat> computeimpliedCovarianceDerivatives(const SEMCpp& SEM,
                                                           const arma::mat& IminusAInverse,
                                                           bool raw){
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
  
  // the following vector will hold the derivatives with respect to each parameter
  std::vector<arma::mat> derivativesOfCovariance(uniqueParameterLabels.size());
  
  // pre-compute some repeatedly used elements
  arma::mat FIminusAInverse = Fmatrix*IminusAInverse;
  arma::mat Aelement1 = Smatrix*arma::trans(IminusAInverse)*arma::trans(Fmatrix);
  arma::mat Aelement2 =  Fmatrix*IminusAInverse*Smatrix;
  
  // iterate over all parameters and compute derivatives
  for(int p = 0; p < uniqueParameterLabels.size(); p++){
    
    if(parameterLocations.at(p).compare("Mvector") == 0) {
      continue; // derivative is zero, as initiated above
    }
    if(parameterLocations.at(p).compare("Smatrix") == 0){

      if(isVariance.at(p) && raw){
        // find value of variance
        double valueIs =-9999.9999;
        for(int j = 0; j < SEM.parameterTable.label.length(); j++){
          std::string compto = Rcpp::as<std::string>(SEM.parameterTable.label.at(j));
          if(uniqueParameterLabels.at(p).compare(compto) == 0){
            valueIs = SEM.parameterTable.value.at(j);
            break;
          }
        }
        if(valueIs == -9999.9999){
          Rcpp::stop("Could not find variance value while computing derivatives.");
        }
        derivativesOfCovariance.at(p) = FIminusAInverse*(valueIs*positionInLocation.at(p))*arma::trans(FIminusAInverse);
      }else{
        derivativesOfCovariance.at(p) = FIminusAInverse*positionInLocation.at(p)*arma::trans(FIminusAInverse);
      }
      
      continue;
      
    }
    if(parameterLocations.at(p).compare("Amatrix") == 0){
      
      derivativesOfCovariance.at(p) = Fmatrix*(
        IminusAInverse*
        positionInLocation.at(p)*
        IminusAInverse)*
        Aelement1 +
        Aelement2*
        arma::trans(IminusAInverse*
        positionInLocation.at(p)*
        IminusAInverse
        )*arma::trans(Fmatrix);
      continue;
    }
    
    Rcpp::stop("Could not locate parameter while computing gradients.");
  }
  return(derivativesOfCovariance);
}


arma::mat computeLogDetSigmas(const SEMCpp& SEM,
                              const arma::mat& IminusAInverse,
                              const std::vector<arma::mat>& derivativesOfCovariance,
                              const std::vector<arma::mat>& impliedCovInverses,
                              bool raw){
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
  
  // the following matrix will hold the derivatives wrt each parameter (in its columns)
  // separated for each missingness pattern in its rows.
  arma::mat logDetSigmas(numberOfMissingnessPatterns, uniqueParameterLabels.size(), arma::fill::zeros);
  
  for(int m = 0; m < numberOfMissingnessPatterns; m++){
    
    for(int p = 0; p < uniqueParameterLabels.size(); p++){
      if(parameterLocations.at(p).compare("Mvector") == 0) {continue;}
      if(parameterLocations.at(p).compare("Smatrix") == 0 || parameterLocations.at(p).compare("Amatrix") == 0){
        
        const arma::mat& currentCovDerivative = derivativesOfCovariance.at(p);
        
        logDetSigmas(m,p) = arma::trace(impliedCovInverses.at(m)*currentCovDerivative(dataSubsets.at(m).notMissing, dataSubsets.at(m).notMissing));
        continue;
        
      }
      Rcpp::stop("Encountered parameter in unknown location");
    }
  }
  
  return(logDetSigmas);
}