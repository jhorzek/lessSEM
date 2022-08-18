#include "SEMWithTransformations.h"

void parametersWithTransformations::asTransformation(Rcpp::StringVector parameterLabels, // labels of parameters which are
                                                     // transformations of other parameters
                                                     SEXP transformationFunctionSEXP)
{
  
  isTransformation.resize(uniqueParameterLabels.length());
  for(int p = 0; p < uniqueParameterLabels.length(); p++)
  { 
    isTransformation.at(p) = false;
  }
  for(int p = 0; p < uniqueParameterLabels.length(); p++)
  {
    for(int t = 0; t < parameterLabels.length(); t++)
    {
      if(uniqueParameterLabels.at(p) == parameterLabels.at(t)) 
      {
        isTransformation.at(p) = true;
        break;
      }
    }
  }
  // create pointer to the transformation function 
  Rcpp::XPtr<transformationFunctionPtr> xpTransformationFunction(transformationFunctionSEXP);
  transformationFunction = *xpTransformationFunction; 
}

void parametersWithTransformations::transform()
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

void SEMWithTransformationsCpp::addTransformation(Rcpp::StringVector parameterLabels, // labels of parameters which are
                                               // transformations of other parameters
                                               SEXP transformationFunctionSEXP)
{
  
  parameterTable.asTransformation(parameterLabels,
                                  transformationFunctionSEXP
  );
}

void SEMWithTransformationsCpp::computeTransformations()
{
  parameterTable.transform();
}

RCPP_MODULE(SEMWithTransformations_cpp){
  using namespace Rcpp;
  Rcpp::class_<SEMWithTransformationsCpp>( "SEMWithTransformationsCpp" )
  // methods
  .method( "addTransformation", &SEMWithTransformationsCpp::addTransformation, "Add a transformation function. Expects parameterLabels and pointer to function.")
  .method( "computeTransformations", &SEMWithTransformationsCpp::computeTransformations, "Compute all transformations")
  ;
}
