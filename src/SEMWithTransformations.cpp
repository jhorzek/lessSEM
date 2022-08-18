#include "SEMWithTransformations.h"

void parametersWithTransformations::asTransformation(SEXP transformationFunctionSEXP)
{
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

void SEMWithTransformationsCpp::addTransformation(SEXP transformationFunctionSEXP)
{
  parameterTable.asTransformation(transformationFunctionSEXP);
}

void SEMWithTransformationsCpp::computeTransformations()
{
  parameterTable.transform();
}

RCPP_MODULE(SEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<SEMCpp>( "SEMCpp" )
    .constructor("Creates a new SEMCpp.")
    .field_readonly( "A", &SEMCpp::Amatrix, "Matrix with directed effects")
    .field_readonly( "S", &SEMCpp::Smatrix, "Matrix with undirected paths")
    .field_readonly( "F", &SEMCpp::Fmatrix, "Filter matrix to separate latent and observed variables")
    .field_readonly( "m", &SEMCpp::Mvector, "Vector with means of observed and latent variables")
    .field_readonly( "impliedCovariance", &SEMCpp::impliedCovariance, "implied covariance matrix")
    .field_readonly( "impliedMeans", &SEMCpp::impliedMeans, "implied means vector")
    .field_readonly( "m2LL", &SEMCpp::m2LL, "minus 2 log-likelihood")
    .field_readonly( "rawData", &SEMCpp::rawData, "raw data set")
    .field_readonly( "manifestNames", &SEMCpp::manifestNames, "names of manifest variables")
    .field_readonly( "wasFit", &SEMCpp::wasFit, "names of manifest variables")
  
  // methods
  .method( "setMatrix", &SEMCpp::setMatrix, "Fills the elements of a model matrix. Expects a char (A, S, or F), and a matrix with values")
  .method( "setVector", &SEMCpp::setVector, "Fills the elements of a model vector. Expects a char (m), and a vector with values")
  .method( "initializeParameters", &SEMCpp::initializeParameters, "Initializes the parameters of the model. Expects StringVector with labels, StringVector with location, uvec with rows, uvec with columns, vec with values, and vec with rawValues")
  .method( "setParameters", &SEMCpp::setParameters, "Changes the parameters of a model. Expects StringVector with labels, vec with values, and boolean indicating if the values are raw (TRUE) or transformed (FALSE).")
  .method( "addDerivativeElement", &SEMCpp::addDerivativeElement, "Add an element to the internal derivative structure. string label, string with location, and boolean indicating if it is a variance, matrix with positions.")
  .method( "addRawData", &SEMCpp::addRawData, "Adds a raw data set. Expects matrix and vector with labels of manifest variables.")
  .method( "addSubset", &SEMCpp::addSubset, "Adds a subset to the data. Expects the sample size N of the subset, the number of observed variables without missings, uvec with indices of notMissing values, matrix with covariance, colvec with means, raw data without missings.")
  .method( "implied", &SEMCpp::implied, "Computes implied means and covariance matrix")
  .method( "fit", &SEMCpp::fit, "Fits the model. Returns -2 log likelihood")
  .method( "getParameters", &SEMCpp::getParameters, "Returns a data frame with model parameters.")
  .method( "getParameterLabels", &SEMCpp::getParameterLabels, "Returns a vector with unique parameter labels as used internally.")
  .method( "getGradients", &SEMCpp::getGradients, "Returns a matrix with scores.")
  .method( "getScores", &SEMCpp::getScores, "Returns a matrix with scores.")
  .method( "getHessian", &SEMCpp::getHessian, "Returns the hessian of the model. Expects the labels of the parameters and the values of the parameters as well as a boolean indicating if these are raw. Finally, a double (eps) controls the precision of the approximation.")
  ;

  Rcpp::class_<SEMWithTransformationsCpp>( "SEMWithTransformationsCpp" )
    .derives<SEMCpp>("SEMCpp")

  // additional methods
  .method( "addTransformation", &SEMWithTransformationsCpp::addTransformation, "Add a transformation function. Expects parameterLabels and pointer to function.")
  .method( "computeTransformations", &SEMWithTransformationsCpp::computeTransformations, "Compute all transformations")
  ;
}
