#include "parameters.h"
#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

void parameters::initialize(Rcpp::StringVector label_,
                            Rcpp::StringVector location_,
                            arma::uvec row_,
                            arma::uvec col_,
                            arma::vec value_,
                            arma::vec rawValue_,
                            std::vector<bool> isTransformation_){
  
  std::string currentLabel;
  // create hash values
  for(int i = 0; i < label_.length(); i++){
    
    currentLabel = Rcpp::as< std::string >(label_.at(i));
    
    if(parameterMap.count(currentLabel) == 0){
      // the parameter does not yet exist
      uniqueParameterLabels.push_back(label_.at(i));
      uniqueParameterValues.push_back(value_.at(i));
      uniqueRawParameterValues.push_back(rawValue_.at(i));
      uniqueParameterLocations.push_back(location_.at(i));
      uniqueIsTransformation.push_back(isTransformation_.at(i));
      
      parameterMap[currentLabel].value = value_.at(i);
      parameterMap[currentLabel].rawValue = rawValue_.at(i);
      parameterMap[currentLabel].changed = true;
      
      parameterMap[currentLabel].location = Rcpp::as< std::string >(location_.at(i));
      
      if((row_.at(i) == col_.at(i)) & (location_.at(i) == "Smatrix")){
        parameterMap[currentLabel].isVariance = true;
      }else{
        parameterMap[currentLabel].isVariance = false;
      }
      
      parameterMap[currentLabel].isTransformation = isTransformation_.at(i);
    }
    
    parameterMap.at(currentLabel).row.push_back(row_.at(i));
    parameterMap.at(currentLabel).col.push_back(col_.at(i));
  }
}

Rcpp::DataFrame parameters::getParameters(){
  Rcpp::LogicalVector isTransformation(uniqueParameterLabels.length());
  for(int i = 0; i < uniqueParameterLabels.length(); i++){
    uniqueParameterValues.at(i) = parameterMap[Rcpp::as< std::string >(uniqueParameterLabels.at(i))].value;
    uniqueRawParameterValues.at(i) = parameterMap[Rcpp::as< std::string >(uniqueParameterLabels.at(i))].rawValue;
    isTransformation.at(i) = parameterMap[Rcpp::as< std::string >(uniqueParameterLabels.at(i))].isTransformation;
  }
  Rcpp::DataFrame parameterFrame =
    Rcpp::DataFrame::create( Rcpp::Named("label") = uniqueParameterLabels, 
                             Rcpp::Named("value") = uniqueParameterValues,
                             Rcpp::Named("rawValue") = uniqueRawParameterValues,
                             Rcpp::Named("location") = uniqueParameterLocations,
                             Rcpp::Named("isTransformation") = isTransformation
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
      if(parameterMap.at(parameterLabel).isTransformation) Rcpp::stop("Cannot change transformed parameters");
      
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

void parameters::addTransformation(SEXP transformationFunctionSEXP, 
                                   Rcpp::List transformationList_)
{
  hasTransformations = true;
  // create pointer to the transformation function 
  Rcpp::XPtr<transformationFunctionPtr> xpTransformationFunction(transformationFunctionSEXP);
  transformationFunction = *xpTransformationFunction; 
  transformationList = transformationList_;
}

void parameters::transform()
{
  Rcpp::NumericVector params(uniqueParameterLabels.length());
  Rcpp::CharacterVector paramLabels(uniqueParameterLabels.length());
  for(int i = 0; i < uniqueParameterLabels.length(); i++){
    params.at(i) = parameterMap[Rcpp::as< std::string >(uniqueParameterLabels.at(i))].rawValue;
    paramLabels.at(i) = uniqueParameterLabels.at(i);
  }
  params.names() = paramLabels;
  
  params = transformationFunction(params, transformationList);
  
  // also change the parameter values in the parameter map; these are the ones that
  // are actually used internally
  std::string parameterLabel;
  for(int p = 0; p < paramLabels.length(); p++){
    parameterLabel = Rcpp::as< std::string >(paramLabels.at(p));
    parameterMap.at(parameterLabel).rawValue = params.at(p);
    if(parameterMap.at(parameterLabel).isVariance){
      parameterMap.at(parameterLabel).value = std::exp(params.at(p));
    }else{
      parameterMap.at(parameterLabel).value = params.at(p);
    }
  }
}

arma::mat parameters::getTransformationGradients(){
  if(!hasTransformations) Rcpp::stop("Does not have transformations.");
  arma::mat currentGradients(nModelParameters, 
                             nRealParameters,
                             arma::fill::zeros);
  Rcpp::NumericVector parameterValues(uniqueParameterLabels.size());
  Rcpp::CharacterVector parameterLabelsRcpp(uniqueParameterLabels.size());
  arma::uvec selectRows(nModelParameters);
  arma::colvec stepForward, stepBackward;
  std::string currentParameter;
  int j = 0;
  for(int i = 0; i < uniqueParameterLabels.size(); i++){
    currentParameter = uniqueParameterLabels.at(i);
    
    parameterValues.at(i) = parameterMap.at(currentParameter).rawValue;
    parameterLabelsRcpp.at(i) = currentParameter;
    if(parameterMap.at(currentParameter).location.compare("transformation") == 0){
      // we want to remove all parameters which are only in transformations and not in the model
      continue;
    }
    selectRows.at(j) = i;
    j++;
  }
  
  parameterValues.names() = parameterLabelsRcpp;
  
  j = 0;
  for(int i = 0; i < parameterValues.size(); i++){
    currentParameter = parameterLabelsRcpp.at(i);
    if(parameterMap.at(currentParameter).isTransformation) continue;
    parameterValues.at(i) += gradientStepSize;
    
    stepForward = Rcpp::as<arma::colvec>(transformationFunction(parameterValues, transformationList)).rows(selectRows);
    
    parameterValues.at(i) -= 2.0*gradientStepSize;
    
    stepBackward = Rcpp::as<arma::colvec>(transformationFunction(parameterValues, transformationList)).rows(selectRows);
    
    parameterValues.at(i) += gradientStepSize;
    
    currentGradients.col(j) = (stepForward - stepBackward)/(2.0*gradientStepSize);
    j++;
  }
  
  return(currentGradients);
}