#ifndef PARAMETERMODULE_H
#define PARAMETERMODULE_H

#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp :: depends ( RcppArmadillo )]]

class parameters{
public: 
  Rcpp::StringVector label; // parameter label
  Rcpp::StringVector location; // name of the matrix where the parameter is located
  arma::uvec row; // row in the matrix 
  arma::uvec col; // col in the matrix
  arma::vec value; // parameter value
  arma::vec rawValue; // raw (untransformed) parameter value
  std::vector<bool> changed; // indicates if a parameter changed its value
  std::map<std::string,std::vector<int>> parameterToLocation; // tells us where 
  // the parameter is located in the vectors above

  // constructor
  parameters(){};
  
  //
  void initialize(Rcpp::StringVector label_,
             Rcpp::StringVector location_,
             arma::uvec row_,
             arma::uvec col_,
             arma::vec value_,
             arma::vec rawValue_);
  
  // getter
  Rcpp::DataFrame getParameters();
  
  // setter
  void setParameters(Rcpp::StringVector label_,
                     arma::vec value_,
                     bool raw
                     );
  
};

#endif