#include "parameters.hpp"
#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

void parameters::initialize(Rcpp::StringVector label_,
                            Rcpp::StringVector location_,
                            arma::uvec row_,
                            arma::uvec col_,
                            arma::vec value_,
                            arma::vec rawValue_){
  
  label = label_;
  location = location_;
  row = row_;
  col = col_;
  value = value_;
  rawValue = rawValue_;
}

Rcpp::DataFrame parameters::getParameters(){
  
  Rcpp::DataFrame parameterFrame =
    Rcpp::DataFrame::create( Rcpp::Named("label") = label, 
                             Rcpp::Named("location") = location,
                             Rcpp::Named("row") = row,
                             Rcpp::Named("col") = col,
                             Rcpp::Named("value") = value,
                             Rcpp::Named("rawValue") = rawValue
    );
  return(parameterFrame);
}

void parameters::setParameters(Rcpp::StringVector label_,
                               arma::vec value_,
                               bool raw){
  bool found = false;
  
  for(int element1 = 0; element1 < label_.size(); element1++){
    // Note: Parameters can occur multiple times; we therefore have
    // to iterate over all labels and can't break early.
    found = false;
    
    for(int element2 = 0; element2 < label.size(); element2++){
      if(label_.at(element1) == label.at(element2)){
        // parameter was found -> change value
        found = true;
        if(!raw){
          // if the parameter values were not provided in raw format,
          // all variances are assumed to be transformations 
          value.at(element2) = value_.at(element1);

          if(row.at(element2) == col.at(element2) && location(element2) == "Smatrix"){

            rawValue.at(element2) = std::log(value_.at(element1));
          }else{
            rawValue.at(element2) = value_.at(element1);
          }
        }else{
          // if the parameter values were provided in raw format,
          // all raw variances are assumed to be unbounded;
          // they will be transformed to get the actual variances
          rawValue.at(element2) = value_.at(element1);

          if(row.at(element2) == col.at(element2) && location(element2) == "Smatrix"){
            value.at(element2) = std::exp(value_.at(element1));
          }else{
            value.at(element2) = value_.at(element1);
          }
        }
      }
    }
    if(!found){
      Rcpp::Rcout << label_.at(element1) << std::endl;
      Rcpp::stop("Parameter not found!"); 
      }
  }
  
  return;
}