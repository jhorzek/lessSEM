#ifndef DERIVATIVEMODULE_H
#define DERIVATIVEMODULE_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

class derivativeElements{
public:
  std::vector<std::string> uniqueLabels;
  std::vector<std::string> uniqueLocations;
  std::vector<bool> isVariance;
  std::vector<arma::mat> positionInLocation; // matrix with 1s in location of the parameter
  std::vector<bool> wasInitialized; // prevent initializing the same parameter twice
  
  // constructor
  derivativeElements(){};
  
  void initialize(int nParam,
                  Rcpp::StringVector uniqueParameterLabels,
                  Rcpp::StringVector uniqueParameterLocations);
    
    void addDerivativeElement(std::string label_, std::string location_, bool isVariance_, arma::mat positionMatrix_);
};

#endif