#ifndef PARAMETERMODULE_H
#define PARAMETERMODULE_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

class parameters{
public: 
  Rcpp::CharacterVector label; // parameter label
  Rcpp::CharacterVector location; // name of the matrix where the parameter is located
  arma::uvec row; // row in the matrix 
  arma::uvec col; // col in the matrix
  arma::vec value; // parameter value
  arma::vec rawValue; // raw (untransformed) parameter value
  
  // constructor
  parameters(){};
  
  //
  void initialize(Rcpp::CharacterVector label_,
             Rcpp::CharacterVector location_,
             arma::uvec row_,
             arma::uvec col_,
             arma::vec value_,
             arma::vec rawValue_);
  
  // getter
  Rcpp::DataFrame getParameters();
  
  // setter
  void setParameters(Rcpp::CharacterVector label_,
                     arma::vec value_,
                     bool raw
                     );
  
};

#endif