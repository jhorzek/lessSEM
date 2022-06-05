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
  changed.resize(label.size());
  std::fill(changed.begin(), changed.end(), true);
  
  // create hashs
  for(int i = 0; i < label.length(); i++){
    
    if(parameterToLocation.find(Rcpp::as< std::string >(label.at(i))) == parameterToLocation.end()){
      std::vector<int> locations;
      locations.push_back(i);
      parameterToLocation[Rcpp::as< std::string >(label.at(i))] = locations;
      continue;
    }else{
      parameterToLocation[Rcpp::as< std::string >(label.at(i))].push_back(i);
    }
  }
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
  
  for(int param = 0; param < label_.size(); param++){
    // all parameters are saved as a hash in parameterToLocation, where
    // the value indicates the locations in the parameter table
    std::vector<int> &locations = parameterToLocation[Rcpp::as< std::string >(label_.at(param))];
    
    // now iterate over these elements and set the parameters
    
    for(int locs = 0; locs < locations.size(); locs++){
      
      if(!raw){
        // if the parameter values were not provided in raw format,
        // all variances are assumed to be transformations 
        if(value.at(locations.at(locs)) == value_.at(param)){
          continue;
        }
        
        value.at(locations.at(locs)) = value_.at(param);
        changed.at(locations.at(locs)) = true;
        
        if(row.at(locations.at(locs)) == col.at(locations.at(locs)) && location(locations.at(locs)) == "Smatrix"){
          
          rawValue.at(locations.at(locs)) = std::log(value_.at(param));
          
        }else{
          
          rawValue.at(locations.at(locs)) = value_.at(param);
          
        }
      }else{
        // if the parameter values were provided in raw format,
        // all raw variances are assumed to be unbounded;
        // they will be transformed to get the actual variances
        
        if(rawValue.at(locations.at(locs)) == value_.at(param)){
          continue;
        }
        
        changed.at(locations.at(locs)) = true;
        
        rawValue.at(locations.at(locs)) = value_.at(param);
        
        if(row.at(locations.at(locs)) == col.at(locations.at(locs)) && location(locations.at(locs)) == "Smatrix"){
          value.at(locations.at(locs)) = std::exp(value_.at(param));
        }else{
          value.at(locations.at(locs)) = value_.at(param);
        }
      }
    }
  }
  
  return;
}