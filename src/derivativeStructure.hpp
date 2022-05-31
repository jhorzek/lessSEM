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
  
  // constructor
  derivativeElements(){};
  
  void addDerivativeElement(std::string label_, std::string location_, bool isVariance_, arma::mat positionMatrix_);
};

#endif