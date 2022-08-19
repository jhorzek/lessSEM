#include <RcppArmadillo.h>
#include "SEM.h"
#include "implied.h"
#include "fit.h"
#include "scores.h"
#include "gradients.h"
#include "hessian.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

//' @name SEMCpp
//' 
//' @title SEMCpp class
//' 
//' @description internal SEM representation
//' 
//' @field new Creates a new SEMCpp.
//' @field A Matrix with directed effects
//' @field S Matrix with undirected paths
//' @field F Filter matrix to separate latent and observed variables
//' @field m Vector with means of observed and latent variables
//' @field impliedCovariance implied covariance matrix
//' @field impliedMeans implied means vector
//' @field m2LL -2 log-likelihood
//' @field rawData raw data set
//' @field manifestNames names of manifest variables
//' @field wasFit names of manifest variables
//' @field setMatrix Fills the elements of a model matrix. Expects a char 
//' (A, S, or F), and a matrix with values
//' @field setVector Fills the elements of a model vector. 
//' Expects a char (m), and a vector with values
//' @field initializeParameters Initializes the parameters of the model. 
//' Expects StringVector with labels, StringVector with location, 
//' uvec with rows, uvec with columns, vec with values, and vec with rawValues
//' @field setParameters Changes the parameters of a model. Expects StringVector 
//' with labels, vec with values, and boolean indicating if the values are raw (TRUE) 
//' or transformed (FALSE).
//' @field addDerivativeElement Add an element to the internal derivative structure. 
//' string label, string with location, and boolean indicating if it is a variance, matrix with positions.
//' @field addRawData Adds a raw data set. Expects matrix and vector with labels of manifest variables.
//' @field addSubset Adds a subset to the data. Expects the sample size N of the subset, 
//' the number of observed variables without missings, uvec with indices of notMissing values, 
//' matrix with covariance, colvec with means, raw data without missings.
//' @field implied Computes implied means and covariance matrix
//' @field fit Fits the model. Returns -2 log likelihood
//' @field getParameters Returns a data frame with model parameters.
//' @field getParameterLabels Returns a vector with unique parameter labels as used internally.
//' @field getGradients Returns a matrix with scores.
//' @field getScores Returns a matrix with scores.
//' @field getHessian Returns the hessian of the model. Expects the labels of the 
//' parameters and the values of the parameters as well as a boolean indicating if 
//' these are raw. Finally, a double (eps) controls the precision of the approximation.


bool SEMCpp::checkModel(){

  int N_ = 0;
  for(int i = 0; i < data.dataSubsets.size(); i++){
    N_ += data.dataSubsets.at(i).N;
  }
  if(N_ != rawData.n_rows){
    Rcpp::stop("The number of subjects in the subsets does not match the rows of the raw data matrix.");
  }
  
  // check if all parameters have defined derivative elements
  std::string currentParameterLabel;
  std::string currentParameterLocation;
  // bool wasFound = false;
  // for(int p = 0; p < parameterTable.label.length(); p++){
  //   currentParameterLabel = Rcpp::as<std::string> (parameterTable.label(p));
  //   currentParameterLocation = Rcpp::as<std::string> (parameterTable.location(p));
  //   wasFound = false;
  //   for(int p2 = 0; p2 < derivElements.uniqueLabels.size(); p2++){
  //     if(currentParameterLabel.compare(derivElements.uniqueLabels.at(p2)) == 0){
  //       if(currentParameterLocation.compare(derivElements.uniqueLocations.at(p2)) != 0){
  //         Rcpp::stop("The location of the parameter in the parameter object and the derivElements does not match.");
  //       }
  //       wasFound = true;
  //     }
  //   }
  // }
  
  return(true);
}

void SEMCpp::addRawData(arma::mat rawData_, Rcpp::StringVector manifestNames_, arma::uvec personInSubset_){
  if((currentStatus != addedDerivatives) & (currentStatus != addedRawData)){
    Rcpp::stop("Please define the model matrices and add parameters as well as add derivative elements before calling addRawData.");
  }
  currentStatus = addedRawData;
  
  rawData = rawData_;
  personInSubset = personInSubset_;
  manifestNames = manifestNames_;
  
  return;
}

void SEMCpp::addSubset(int N_,
                       arma::uvec persons_, // indices for which persons are in this subset
                       int observed_, // number of variables without missings
                       arma::uvec notMissing_, // vector with indices of non-missing values
                       // the following elements are only relevant for N>1
                       arma::mat covariance_,
                       arma::colvec means_,
                       // raw data is required for N == 1
                       arma::mat rawData_){
  
  if((currentStatus != addedRawData) & (currentStatus != addedSubsets)){
    Rcpp::stop("Please define the model matrices and add parameters as well as add derivative elements before calling addRawData.");
  }
  currentStatus = addedSubsets;
  
  data.addSubset(N_,
                 persons_,
                 observed_, // number of variables without missings
                 notMissing_, // vector with indices of non-missing values
                 // the following elements are only relevant for N>1
                 covariance_,
                 means_,
                 // raw data is required for N == 1
                 rawData_);
  
  return;
  
}

void SEMCpp::removeSubset(int whichSubset){
  data.removeSubset(whichSubset);
  return;
}

void SEMCpp::setMatrix(std::string whichMatrix, arma::mat values){
  
  currentStatus = addedMatrices; 
  
  if(whichMatrix.compare("A") == 0){
    Amatrix = values;
    return;
  }
  if(whichMatrix.compare("S") == 0){
    Smatrix = values;
    return;
  }
  if(whichMatrix.compare("F") == 0){
    Fmatrix = values;
    return;
  }
  Rcpp::stop("Unknown whichMatrix. Can only assign A, S, or F.");
  
}

void SEMCpp::setVector(std::string whichVector, arma::colvec values){
  
  if(whichVector.compare("m") == 0){
    Mvector = values;
    return;
  }
  
  Rcpp::stop("Unknown whichVector. Can only assign m.");
}

void SEMCpp::initializeParameters(Rcpp::StringVector label_,
                                  Rcpp::StringVector location_,
                                  arma::uvec row_,
                                  arma::uvec col_,
                                  arma::vec value_,
                                  arma::vec rawValue_,
                                  std::vector<bool> isTransformation){
  if(currentStatus != addedMatrices){
    Rcpp::stop("Please define the model matrices before adding parameters");
  }
  currentStatus = addedParameters;
  
  parameterTable.initialize(label_,
                            location_,
                            row_,
                            col_,
                            value_,
                            rawValue_,
                            isTransformation);
  
  // also initialize the derivative elements
  parameterTable.nModelParameters = 0;
  parameterTable.nTransformationParameters = 0;
  std::string currentParameter;
  for(int i = 0; i < parameterTable.uniqueParameterLabels.length(); i++){
    currentParameter = parameterTable.uniqueParameterLabels.at(i);
    if(parameterTable.parameterMap.at(currentParameter).location.compare("transformation") == 0){
      // parameter is not in the model but only in the transformations
      parameterTable.nTransformationParameters++;
      continue;
    }
    if(parameterTable.parameterMap.at(currentParameter).isTransformation){
      // parameter is in the model, but is a function of other parameters
      parameterTable.nModelParameters++;
      continue;
    }
    parameterTable.nTransformationParameters++;
    parameterTable.nModelParameters++;
  }

  derivElements.initialize(
    parameterTable.nModelParameters,
    parameterTable.uniqueParameterLabels,
    parameterTable.uniqueParameterLocations
  );
  
  return;
}

void SEMCpp::addDerivativeElement(std::string label_, 
                                  std::string location_, 
                                  bool isVariance_, 
                                  arma::mat positionMatrix_){
  if(currentStatus != addedParameters & currentStatus != addedDerivatives){
    Rcpp::stop("Please define the model matrices and add parameters before calling addDerivativeElement.");
  }
  currentStatus = addedDerivatives;
  
  derivElements.addDerivativeElement(label_, 
                                     location_, 
                                     isVariance_, 
                                     positionMatrix_);
  return;
}

void SEMCpp::addTransformation(SEXP transformationFunctionSEXP)
{
  hasTransformations = true;
  parameterTable.addTransformation(transformationFunctionSEXP);
}

void SEMCpp::computeTransformations()
{
  parameterTable.transform();
}

void SEMCpp::setParameters(Rcpp::StringVector label_,
                           arma::vec value_,
                           bool raw){
  currentStatus = changedParameters;
  
  if(parameterTable.hasTransformations & (!raw)) 
    Rcpp::stop("SEMs with transformation currently only support setParameters with raw = true.");
  
  wasFit = false; // reset fit 
  // step one: change parameters in parameterTable
  parameterTable.setParameters(label_, value_, raw);
  if(parameterTable.hasTransformations)
   parameterTable.transform();
  // step two: change parameters in model matrices
  
  for (auto const& param : parameterTable.parameterMap)
  {
    if(!param.second.changed) continue;
    
    if(param.second.location.compare("Amatrix") == 0){
      parameterTable.AChanged = true;
      
      for(int elem = 0; elem < param.second.row.size(); elem ++){
        
        Amatrix(param.second.row.at(elem), 
                param.second.col.at(elem)) = param.second.value;
      }
      continue;
    }
    
    if(param.second.location.compare("Smatrix") == 0){
      parameterTable.SChanged = true;
      
      for(int elem = 0; elem < param.second.row.size(); elem ++){
        
        Smatrix(param.second.row.at(elem), 
                param.second.col.at(elem)) = param.second.value;
      }
      continue;
    }
    
    if(param.second.location.compare("Mvector") == 0){
      parameterTable.mChanged = true;
      
      for(int elem = 0; elem < param.second.row.size(); elem ++){
        
        Mvector(param.second.row.at(elem), 
                param.second.col.at(elem)) = param.second.value;
        
      }
      continue;
    }
    
    if(param.second.location.compare("transformation") == 0){
      continue;
    }
    
    Rcpp::stop("NOT FOUND");
  }
  
  // reset all fit elements
  m2LL = 0.0;
  if(parameterTable.AChanged | parameterTable.SChanged)  impliedCovariance.fill(arma::datum::nan);
  if(parameterTable.mChanged) impliedMeans.fill(arma::datum::nan);
  return;
}

// getter
Rcpp::DataFrame SEMCpp::getParameters(){
  
  return(parameterTable.getParameters());
  
}

Rcpp::StringVector SEMCpp::getParameterLabels(){
  
  return(parameterTable.uniqueParameterLabels);
  
}
// fit functions


// [[Rcpp :: depends ( RcppArmadillo )]]

void SEMCpp::implied(){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  currentStatus = computedImplied;

  // step one: compute implied mean and covariance
  if(parameterTable.AChanged | parameterTable.SChanged) impliedCovariance = computeImpliedCovariance(Fmatrix, Amatrix, Smatrix);
  if(parameterTable.mChanged | parameterTable.AChanged) impliedMeans = computeImpliedMeans(Fmatrix, Amatrix, Mvector);
  
  parameterTable.AChanged = false;
  parameterTable.SChanged = false;
  parameterTable.mChanged = false;
  
  return;
}


double SEMCpp::fit(){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  currentStatus = fitted;
  
  m2LL = 0.0;
  wasFit = true;
  
  // step one: compute implied mean and covariance
  implied();
  
  // step two: compute fit for each subset
  
  for(int s = 0; s < data.nGroups; s++){
    
    // convenient pointer to current subset:
    subset& currentSubset = data.dataSubsets.at(s);
    arma::mat currentImpliedCovariance = impliedCovariance.submat(currentSubset.notMissing, currentSubset.notMissing);
    arma::colvec currentImpliedMeans = impliedMeans.rows(currentSubset.notMissing);
    
    if(currentSubset.N == 1){
      
      currentSubset.m2LL = computeIndividualM2LL(currentSubset.observed, 
                                                 currentSubset.rawData(currentSubset.notMissing), 
                                                 currentImpliedMeans, 
                                                 currentImpliedCovariance);
      m2LL += currentSubset.m2LL;
      
    }else if(currentSubset.N > 1){
      
      currentSubset.m2LL = computeGroupM2LL(currentSubset.N, 
                                            currentSubset.observed, 
                                            currentSubset.means, 
                                            currentSubset.covariance,
                                            currentImpliedMeans, 
                                            currentImpliedCovariance);
      m2LL += currentSubset.m2LL;
      
    }else{
      Rcpp::stop("Unknown N in data subset");
    }
    
  }
  
  return(m2LL);
}

arma::mat SEMCpp::getScores(bool raw){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  if((currentStatus != computedImplied) & (currentStatus != fitted)){
    Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  }
  if(hasTransformations) Rcpp::stop("Not yet implemented for models with transformations.");
  arma::mat scoresMat = scores(*this, raw);
  
  return(scoresMat);
}

arma::rowvec SEMCpp::getGradients(bool raw){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  if((currentStatus != computedImplied) & (currentStatus != fitted)){
    Rcpp::stop("The model implied matrices have not been computed yet. Call Model$implied() first.");
  }
  arma::rowvec gradients = gradientsByGroup(*this, raw);
  
  if(hasTransformations){
    if(!raw) Rcpp::stop("Gradients with raw = false currently not supported when using transformations.");
    // in case of transformations, we can make use of the chain rule to get the derivative
    // with respect to the true underlying parameters. gradientsByGroup will only return
    // the gradients of the transformed parameters in the SEM 
    arma::mat transformationGradients = parameterTable.getTransformationGradients();
    return(gradients*transformationGradients);
  }
  
  return(gradients);
}

arma::mat SEMCpp::getHessian(Rcpp::StringVector label_,
                             arma::vec value_,
                             bool raw,
                             double eps){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  if((currentStatus != computedImplied) & (currentStatus != fitted)){
    Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  }
  
  if(hasTransformations){
    Rcpp::stop("Hessian for transformed parameters is not yet implemented");
  }
  
  arma::mat hessian = approximateHessian(*this, 
                                         label_,
                                         value_,
                                         raw,
                                         eps);
  
  return(hessian);
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
  .method( "addTransformation", &SEMCpp::addTransformation, "Add a transformation function. Expects parameterLabels and pointer to function.")
  .method( "computeTransformations", &SEMCpp::computeTransformations, "Compute all transformations")
  ;
}
