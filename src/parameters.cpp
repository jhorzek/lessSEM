#include "parameters.h"
#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

void parameters::initialize(Rcpp::StringVector label_,
                            Rcpp::StringVector location_,
                            arma::uvec row_,
                            arma::uvec col_,
                            arma::vec value_,
                            arma::vec rawValue_){
  
  std::string currentLabel;
  // create hash values
  for(int i = 0; i < label_.length(); i++){
    
    currentLabel = Rcpp::as< std::string >(label_.at(i));
    
    if(parameterMap.count(currentLabel) == 0){
      // the parameter does not yet exist
      uniqueParameterLabels.push_back(label_.at(i));
      uniqueParameterValues.push_back(value_.at(i));
      uniqueRawParameterValues.push_back(rawValue_.at(i));
      
      parameterMap[currentLabel].value = value_.at(i);
      parameterMap[currentLabel].rawValue = rawValue_.at(i);
      parameterMap[currentLabel].changed = true;
      
      parameterMap[currentLabel].location = Rcpp::as< std::string >(location_.at(i));
      
      if((row_.at(i) == col_.at(i)) & (location_.at(i) == "Smatrix")){
        parameterMap[currentLabel].isVariance = true;
      }else{
        parameterMap[currentLabel].isVariance = false;
      }
      
    }
    
    parameterMap.at(currentLabel).row.push_back(row_.at(i));
    parameterMap.at(currentLabel).col.push_back(col_.at(i));
  }
}

Rcpp::DataFrame parameters::getParameters(){
  for(int i = 0; i < uniqueParameterLabels.length(); i++){
    uniqueParameterValues.at(i) = parameterMap[Rcpp::as< std::string >(uniqueParameterLabels.at(i))].value;
    uniqueRawParameterValues.at(i) = parameterMap[Rcpp::as< std::string >(uniqueParameterLabels.at(i))].rawValue;
  }
  Rcpp::DataFrame parameterFrame =
    Rcpp::DataFrame::create( Rcpp::Named("label") = uniqueParameterLabels, 
                             Rcpp::Named("value") = uniqueParameterValues,
                             Rcpp::Named("rawValue") = uniqueRawParameterValues
    );
  return(parameterFrame);
}

void parameters::setParameters(Rcpp::StringVector label_,
                               arma::vec value_,
                               bool raw){
  std::string parameterLabel;
  for(int param = 0; param < label_.size(); param++){
    parameterLabel = Rcpp::as< std::string >(label_.at(param));
    
    // all parameters are saved as a hash in parameterMap
    if(!raw){
      if(parameterMap.at(parameterLabel).value == value_.at(param)) continue;
      
      parameterMap.at(parameterLabel).changed = true;
      
      parameterMap.at(parameterLabel).value = value_.at(param);
      
      if(parameterMap.at(parameterLabel).isVariance){
        
        parameterMap.at(parameterLabel).rawValue = std::log(value_.at(param));
        
      }else{
        
        parameterMap.at(parameterLabel).rawValue = value_.at(param);
        
      }
    }else{
      
      if(parameterMap.at(parameterLabel).rawValue == value_.at(param)) continue;
      
      parameterMap.at(parameterLabel).changed = true;
      
      parameterMap.at(parameterLabel).rawValue = value_.at(param);
      
      if(parameterMap.at(parameterLabel).isVariance){
        
        parameterMap.at(parameterLabel).value = std::exp(value_.at(param));
        
      }else{
        
        parameterMap.at(parameterLabel).value = value_.at(param);
        
      }
    }
    
    if(parameterMap.at(parameterLabel).location.compare("Smatrix")) SChanged = true;
    if(parameterMap.at(parameterLabel).location.compare("Amatrix")) AChanged = true;
    if(parameterMap.at(parameterLabel).location.compare("mVector")) mChanged = true;
  }
  
  return;
}

void parameters::asTransformation(SEXP transformationFunctionSEXP)
{
  // create pointer to the transformation function 
  Rcpp::XPtr<transformationFunctionPtr> xpTransformationFunction(transformationFunctionSEXP);
  transformationFunction = *xpTransformationFunction; 
}

void parameters::transform()
{
  for(int i = 0; i < uniqueParameterLabels.length(); i++){
    uniqueRawParameterValues.at(i) = parameterMap[Rcpp::as< std::string >(uniqueParameterLabels.at(i))].rawValue;
  }
  uniqueRawParameterValues.names() = uniqueParameterLabels;
  
  uniqueRawParameterValues = transformationFunction(uniqueRawParameterValues);
  
  // also change the parameter values in the parameter map; these are the ones that
  // are actually used internally
  std::string parameterLabel;
  for(int p = 0; p < uniqueParameterLabels.length(); p++){
    parameterLabel = Rcpp::as< std::string >(uniqueParameterLabels.at(p));
    parameterMap.at(parameterLabel).rawValue = uniqueRawParameterValues.at(p);
    if(parameterMap.at(parameterLabel).isVariance){
      parameterMap.at(parameterLabel).value = std::log(uniqueRawParameterValues.at(p));
    }else{
      parameterMap.at(parameterLabel).value = uniqueRawParameterValues.at(p);
    }
  }
}