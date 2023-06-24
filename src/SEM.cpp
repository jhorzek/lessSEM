#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "SEM.h"
#include "implied.h"
#include "impliedDerivatives.h"
#include "multivariateNormal.h"
#include "WLS.h"
#include "scores.h"
#include "gradients.h"
#include "hessian.h"

// [[Rcpp :: depends ( RcppArmadillo , RcppParallel)]]

//' @name SEMCpp
//' 
//' @title SEMCpp class
//' 
//' @description internal SEM representation
//' 
//' @field new Creates a new SEMCpp.
//' @field fill fills the SEM with the elements from an Rcpp::List
//' @field addTransformation adds transforamtions to a model
//' @field implied Computes implied means and covariance matrix
//' @field fit Fits the model. Returns objective value of the fitting function
//' @field getParameters Returns a data frame with model parameters.
//' @field getEstimator returns the estimator used in the model (e.g., fiml)
//' @field getParameterLabels Returns a vector with unique parameter labels as used internally.
//' @field getGradients Returns a matrix with scores.
//' @field getScores Returns a matrix with scores.
//' @field getHessian Returns the hessian of the model. Expects the labels of the 
//' parameters and the values of the parameters as well as a boolean indicating if 
//' these are raw. Finally, a double (eps) controls the precision of the approximation.
//' @field computeTransformations compute the transformations.
//' @field setTransformationGradientStepSize change the step size of the gradient computation for the transformations
 
 void SEMCpp::fill(Rcpp::List SEMList){
   
   currentStatus = initialized;
   
   // Estimator
   std::string estim_ = SEMList["estimator"];
   if(estim_.compare("fiml") == 0){
     estim = fiml;
   }else if(estim_.compare("wls") == 0){
     estim = wls;
     
     WLSElements.weightsInverse = Rcpp::as<arma::mat>(SEMList["WLSWeightsInverse"]);
     WLSElements.meanstructure = Rcpp::as<bool>(SEMList["WLSMeanstructure"]);
     
   }else{
     Rcpp::stop("Unknown estimator selected. Possible are FIML and WLS.");
   }
   
   // Add all matrices with matrices
   Rcpp::List matrices = SEMList["matrices"];
   
   Amatrix = Rcpp::as<arma::mat>(matrices["A"]);
   Smatrix = Rcpp::as<arma::mat>(matrices["S"]);
   Fmatrix = Rcpp::as<arma::mat>(matrices["F"]);
   Mvector = Rcpp::as<arma::colvec>(matrices["M"]);
   
   // add parameters to the model
   Rcpp::List param = SEMList["parameters"];
   
   Rcpp::StringVector label = param["label"];
   Rcpp::StringVector location = param["location"];
   arma::uvec row = param["row"];
   arma::uvec col = param["col"];
   arma::vec value = param["value"];
   arma::vec rawValue = param["rawValue"];
   std::vector<bool> isTransformation = param["isTransformation"];
   parameterTable.initialize(label,
                             location,
                             row,
                             col,
                             value,
                             rawValue,
                             isTransformation);
   
   // also _initialize_ the derivative elements
   parameterTable.nModelParameters = 0;
   parameterTable.nRealParameters = 0;
   std::string currentParameter;
   for(unsigned int i = 0; i < parameterTable.uniqueParameterLabels.length(); i++){
     currentParameter = parameterTable.uniqueParameterLabels.at(i);
     if(parameterTable.parameterMap.at(currentParameter).location.compare("transformation") == 0){
       // parameter is not in the model but only in the transformations
       parameterTable.nRealParameters++;
       continue;
     }
     if(parameterTable.parameterMap.at(currentParameter).isTransformation){
       // parameter is in the model, but is a function of other parameters
       parameterTable.nModelParameters++;
       continue;
     }
     parameterTable.nRealParameters++;
     parameterTable.nModelParameters++;
   }
   
   derivElements.initialize(
     parameterTable.nModelParameters,
     parameterTable.uniqueParameterLabels,
     parameterTable.uniqueParameterLocations
   );
   
   // after initializing the derivative elements, we now have to fill them
   Rcpp::List fillDerivElementsWith = SEMList["DerivativeElements"];
   
   for(unsigned int de = 0; de < fillDerivElementsWith.length(); de++){
     
     Rcpp::List currentDerivElement = fillDerivElementsWith[de];
     
     std::string label_ = currentDerivElement["label"]; 
     std::string location_ = currentDerivElement["location"];
     bool isVariance_ = currentDerivElement["isVariance"]; 
     arma::mat positionMatrix_ = currentDerivElement["positionMatrix"];
     
     derivElements.addDerivativeElement(label_, 
                                        location_, 
                                        isVariance_, 
                                        positionMatrix_);
   }
   
   // add raw data
   Rcpp::List dataset = SEMList["rawData"];
   rawData = Rcpp::as<arma::mat>(dataset["rawData"]);
   personInSubset = Rcpp::as<arma::uvec>(dataset["personInSubset"]);
   manifestNames = dataset["manifestNames"];
   sampleSize = rawData.n_rows;
   
   // add subsets
   Rcpp::List subsets = SEMList["subsets"];
   
   for(unsigned int s = 0; s < subsets.length(); s++){
     
     Rcpp::List subset = subsets[s];
     
     int N_ = subset["N"];
     arma::uvec persons_  = subset["persons"]; // indices for which persons are in this subset
     int observed_ = subset["observed"]; // number of variables without missings
     arma::uvec notMissing_ = subset["notMissing"]; // vector with indices of non-missing values
     // the following elements are only relevant for N>1
     arma::mat covariance_ = subset["covariance"];
     arma::colvec means_ = subset["means"];
     // raw data is required for N == 1
     arma::mat rawData_ = subset["rawData"];
     
     data.addSubset(N_,
                    persons_,
                    observed_, // number of variables without missings
                    notMissing_, // vector with indices of non-missing values
                    // the following elements are only relevant for N>1
                    covariance_,
                    means_,
                    // raw data is required for N == 1
                    rawData_);
   }
   // intialize some elements for the subsets:
   subsetImpliedMeans.resize(subsets.length());
   subsetImpliedCovariance.resize(subsets.length());
   subsetImpliedCovarianceInverse.resize(subsets.length());
   
   return;
   
 }

void SEMCpp::addTransformation(SEXP transformationFunctionSEXP, 
                               Rcpp::List transformationList)
{
  hasTransformations = true;
  parameterTable.addTransformation(transformationFunctionSEXP, transformationList);
}

bool SEMCpp::checkModel(){
  unsigned int N_ = 0;
  for(unsigned int i = 0; i < data.dataSubsets.size(); i++){
    N_ += data.dataSubsets.at(i).N;
  }
  
  if(N_ != rawData.n_rows){
    Rcpp::stop("The number of subjects in the subsets does not match the rows of the raw data matrix.");
  }
  
  return(true);
  
}


void SEMCpp::computeTransformations()
{
  parameterTable.transform();
}

void SEMCpp::setParameters(Rcpp::StringVector label_,
                           arma::vec value_,
                           bool raw){
  currentStatus = changedParameters;
  
  // if(parameterTable.hasTransformations & (!raw)) 
  //   Rcpp::stop("SEMs with transformation currently only support setParameters with raw = true.");
  
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
      
      for(unsigned int elem = 0; elem < param.second.row.size(); elem ++){
        
        Amatrix(param.second.row.at(elem), 
                param.second.col.at(elem)) = param.second.value;
      }
      continue;
    }
    
    if(param.second.location.compare("Smatrix") == 0){
      parameterTable.SChanged = true;
      
      for(unsigned int elem = 0; elem < param.second.row.size(); elem ++){
        
        Smatrix(param.second.row.at(elem), 
                param.second.col.at(elem)) = param.second.value;
      }
      continue;
    }
    
    if(param.second.location.compare("Mvector") == 0){
      parameterTable.mChanged = true;
      
      for(unsigned int elem = 0; elem < param.second.row.size(); elem ++){
        
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
  objectiveValue = 0.0;
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

void SEMCpp::implied(){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  currentStatus = computedImplied;
  
  // step one: compute implied mean and covariance
  if(parameterTable.AChanged){
    IminusAInverse = arma::inv(arma::eye(arma::size(Amatrix))- Amatrix);
  }
  if(parameterTable.AChanged | parameterTable.SChanged) {
    impliedCovarianceFull = computeImpliedCovarianceFull(Amatrix, Smatrix, IminusAInverse);
    impliedCovariance = computeImpliedCovariance(Fmatrix, impliedCovarianceFull);
  }
  if(parameterTable.mChanged | parameterTable.AChanged){
    impliedMeansFull = computeImpliedMeansFull(Amatrix, Mvector, IminusAInverse);
    impliedMeans = computeImpliedMeans(Fmatrix, impliedMeansFull);
  } 
  
  parameterTable.AChanged = false;
  parameterTable.SChanged = false;
  parameterTable.mChanged = false;
  
  return;
}

bool SEMCpp::impliedIsPD(){
  return(impliedCovariance.is_sympd());
}

double SEMCpp::fit(){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  currentStatus = fitted;
  
  objectiveValue = 0.0;
  wasFit = true;
  functionCalls++;
  
  // step one: compute implied mean and covariance
  implied();
  
  double logDetImplied;
  
  // step two: compute fit for each subset
  
  for(unsigned int s = 0; s < data.nGroups; s++){
    
    // convenient pointer to current subset:
    subset& currentSubset = data.dataSubsets.at(s);
    
    subsetImpliedMeans.at(s) = impliedMeans.rows(currentSubset.notMissing);
    subsetImpliedCovariance.at(s) = impliedCovariance.submat(currentSubset.notMissing, currentSubset.notMissing);
    
    if(estim == fiml){
      
      subsetImpliedCovarianceInverse.at(s) =  arma::inv(subsetImpliedCovariance.at(s));
      
      logDetImplied = std::log(arma::det(subsetImpliedCovariance.at(s)));
      
      if(currentSubset.N == 1){
        
        currentSubset.objectiveValue = m2LLMultiVariateNormal(
          currentSubset.rawData(currentSubset.notMissing), 
          subsetImpliedMeans.at(s),
          subsetImpliedCovariance.at(s),
          subsetImpliedCovarianceInverse.at(s),
          logDetImplied);
        objectiveValue += currentSubset.objectiveValue;
        
      }else if(currentSubset.N > 1){
        
        currentSubset.objectiveValue = m2LLGroupMultiVariateNormal(
          currentSubset.N, 
          currentSubset.means,
          subsetImpliedMeans.at(s), 
          currentSubset.covariance,
          subsetImpliedCovariance.at(s),
          subsetImpliedCovarianceInverse.at(s),
          logDetImplied);
        objectiveValue += currentSubset.objectiveValue;
        
      }else{
        Rcpp::stop("Unknown N in data subset");
      }
    }else if (estim == wls){
      
      if(!WLSElements.meanstructure){
        currentSubset.objectiveValue = WLS(
          WLSElements.weightsInverse,
          currentSubset.covariance,
          subsetImpliedCovariance.at(s));
        
      }else{
        currentSubset.objectiveValue = WLS(
          WLSElements.weightsInverse,
          currentSubset.means,
          subsetImpliedMeans.at(s), 
          currentSubset.covariance,
          subsetImpliedCovariance.at(s));
        
      }
      
      // to stay consistent with lavaan, we have to use the following
      // scaling:
      //currentSubset.objectiveValue *=  0.5 * (currentSubset.N-1.0)/currentSubset.N;
      
      // We scale with N because the penalty is also scaled with N.
      currentSubset.objectiveValue *=  currentSubset.N;
      
      objectiveValue += currentSubset.objectiveValue;
      
    }else{
      Rcpp::stop("Unknown estimator");
    }
  }
  
  return(objectiveValue);
}

arma::mat SEMCpp::getScores(bool raw){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  if((currentStatus != computedImplied) && (currentStatus != fitted)){
    Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  }
  if(hasTransformations) Rcpp::stop("Not yet implemented for models with transformations.");
  if(estim != fiml) Rcpp::stop("Currently only implemented for fiml estimation");
  
  if(!detivativesInitialized){
    // there are some elements that are used repeatedly that 
    // are only initialized if the gradients have been computed at least once.
    getGradients(raw);
  }
  
  arma::mat scoresMat = scores(*this, raw);
  
  return(scoresMat);
}

arma::rowvec SEMCpp::getGradients(bool raw){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  if((currentStatus != computedImplied) && (currentStatus != fitted)){
    Rcpp::stop("The model implied matrices have not been computed yet. Call Model$fit() first.");
  }
  gradientCalls++;
  
  // initialize some elements that are used when computing the gradients:
  initializeGradients initializedGradients(*this, raw);
  
  RcppParallel::parallelFor(0, 
                            derivElements.uniqueLabels.size(), //number of parameters
                            initializedGradients
  );
  
  // compute the actual gradients
  gradients = gradientsByGroup(*this, raw);
  
  if(hasTransformations){
    if(!raw) Rcpp::stop("Gradients with raw = false currently not supported when using transformations.");
    // in case of transformations, we can make use of the chain rule to get the derivative
    // with respect to the true underlying parameters. gradientsByGroup will only return
    // the gradients of the transformed parameters in the SEM 
    transformationGradients = parameterTable.getTransformationGradients();
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
  if((currentStatus != computedImplied) && (currentStatus != fitted)){
    Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  }
  
  arma::mat hessian = approximateHessian(*this, 
                                         label_,
                                         value_,
                                         raw,
                                         eps);
  
  return(hessian);
}

void SEMCpp::setTransformationGradientStepSize(double gradientStepSize){
  parameterTable.gradientStepSize = gradientStepSize;
}

std::string SEMCpp::getEstimator(){
  switch(estim){
  case fiml:
    return("fiml");
  case wls:
    return("wls");
  default:
    Rcpp::stop("Cannot find estimator");
  }
}

RCPP_MODULE(SEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<SEMCpp>( "SEMCpp" )
    .constructor("Creates a new SEMCpp.")
  // fields
  .field_readonly("sampleSize", &SEMCpp::sampleSize, "N")
  .field_readonly( "A", &SEMCpp::Amatrix, "Matrix with directed effects")
  .field_readonly( "S", &SEMCpp::Smatrix, "Matrix with undirected paths")
  .field_readonly( "F", &SEMCpp::Fmatrix, "Filter matrix to separate latent and observed variables")
  .field_readonly( "m", &SEMCpp::Mvector, "Vector with means of observed and latent variables")
  .field_readonly( "impliedCovariance", &SEMCpp::impliedCovariance, "implied covariance matrix")
  .field_readonly( "impliedMeans", &SEMCpp::impliedMeans, "implied means vector")
  .field_readonly( "objectiveValue", &SEMCpp::objectiveValue, "minus 2 log-likelihood")
  .field_readonly( "rawData", &SEMCpp::rawData, "raw data set")
  .field_readonly( "manifestNames", &SEMCpp::manifestNames, "names of manifest variables")
  .field_readonly( "wasFit", &SEMCpp::wasFit, "names of manifest variables")
  .field_readonly( "functionCalls", &SEMCpp::functionCalls, "Counts how often the fit function was called")
  .field_readonly( "gradientCalls", &SEMCpp::gradientCalls, "Counts how often the gradient function was called")
  
  // methods
  .method( "fill", &SEMCpp::fill, "Fill all elements of the SEM. Expects and Rcpp::List")
  .method( "implied", &SEMCpp::implied, "Computes implied means and covariance matrix")
  .method( "fit", &SEMCpp::fit, "Fits the model. Returns objective value of the fitting function")
  .method( "setParameters", &SEMCpp::setParameters, "Set the parameters of a model.")
  .method( "getEstimator", &SEMCpp::getEstimator, "returns the used estimator as string.")
  .method( "getParameters", &SEMCpp::getParameters, "Returns a data frame with model parameters.")
  .method( "getParameterLabels", &SEMCpp::getParameterLabels, "Returns a vector with unique parameter labels as used internally.")
  .method( "getGradients", &SEMCpp::getGradients, "Returns a matrix with scores.")
  .method( "getScores", &SEMCpp::getScores, "Returns a matrix with scores.")
  .method( "getHessian", &SEMCpp::getHessian, "Returns the hessian of the model. Expects the labels of the parameters and the values of the parameters as well as a boolean indicating if these are raw. Finally, a double (eps) controls the precision of the approximation.")
  .method( "addTransformation", &SEMCpp::addTransformation, "Add a transformation function. Expects parameterLabels and pointer to function.")
  .method( "computeTransformations", &SEMCpp::computeTransformations, "Compute all transformations")
  .method( "setTransformationGradientStepSize", &SEMCpp::setTransformationGradientStepSize, "Change the step size used in the computation of the transformation gradients.")
  ;
}
