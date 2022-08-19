#include "derivativeStructure.h"
#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

void derivativeElements::initialize(Rcpp::StringVector uniqueParameterLabels,
                                    Rcpp::StringVector uniqueParameterLocations)
  {
  
  uniqueLabels.resize(uniqueParameterLabels.length());
  uniqueLocations.resize(uniqueParameterLabels.length());
  wasInitialized.resize(uniqueParameterLabels.length());
  isVariance.resize(uniqueParameterLabels.length());
  positionInLocation.resize(uniqueParameterLabels.length());
  
  for(int i = 0; i < uniqueParameterLabels.length(); i++){
    uniqueLabels.at(i) = uniqueParameterLabels.at(i);
    uniqueLocations.at(i) = uniqueParameterLocations.at(i);
    wasInitialized.at(i) = false;
  }
  
  return;
}

void derivativeElements::addDerivativeElement(std::string label_, 
                                              std::string location_, 
                                              bool isVariance_, 
                                              arma::mat positionMatrix_){
  
  for(int i = 0; i < uniqueLabels.size(); i++){
    if((uniqueLabels.at(i).compare(label_) == 0) & wasInitialized.at(i)){
      Rcpp::stop("The label passed to addDerivativeElement is already present in the vector uniqueLabels");
    }else if((uniqueLabels.at(i).compare(label_) == 0)){
      wasInitialized.at(i) = true;
      
      if(uniqueLocations.at(i).compare(location_) != 0) {
        Rcpp::stop("Reinitialization with different location.");
      }
      isVariance.at(i) = isVariance_;
      positionInLocation.at(i) = positionMatrix_; // matrix with 1s in location of the parameter
      
      return;
    }else{
      continue;
    }
  }
  
  Rcpp::stop("Could not find parameter in specified parameter labels.");
}