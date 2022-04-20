#include "derivativeStructure.hpp"
#include <RcppArmadillo.h>
#include "config.hpp"
// [[Rcpp :: depends ( RcppArmadillo )]]

void derivativeElements::addDerivativeElement(std::string label_, 
                                              std::string location_, 
                                              bool isVariance_, 
                                              arma::mat positionMatrix_){
  
  for(int i = 0; i < uniqueLabels.size(); i++){
    if(uniqueLabels.at(i).compare(label_) == 0){
      Rcpp::stop("The label passed to addDerivativeElement is already present in the vector uniqueLabels");
    }
  }
  
  uniqueLabels.push_back(label_);
  uniqueLocations.push_back(location_);
  isVariance.push_back(isVariance_);
  positionInLocation.push_back(positionMatrix_); // matrix with 1s in location of the parameter
  
  return;
}