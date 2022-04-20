#include <RcppArmadillo.h>
#include "config.hpp"
#include "SEM.hpp"
#include "implied.hpp"
#include "fit.hpp"
#include "scores.hpp"
#include "gradients.hpp"
#include "hessian.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

bool SEMCpp::checkModel(){
  // some rudimentary model checks
  if(currentStatus < 5){
    Rcpp::stop("The model is incomplete.");
  }
  
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
  bool wasFound = false;
  for(int p = 0; p < parameterTable.label.length(); p++){
    currentParameterLabel = Rcpp::as<std::string> (parameterTable.label(p));
    currentParameterLocation = Rcpp::as<std::string> (parameterTable.location(p));
    wasFound = false;
    for(int p2 = 0; p2 < derivElements.uniqueLabels.size(); p2++){
      if(currentParameterLabel.compare(derivElements.uniqueLabels.at(p2)) == 0){
        if(currentParameterLocation.compare(derivElements.uniqueLocations.at(p2)) != 0){
          Rcpp::stop("The location of the parameter in the parameter object and the derivElements does not match.");
        }
        wasFound = true;
      }
    }
  }
  
  return(true);
}

void SEMCpp::addRawData(arma::mat rawData_, Rcpp::StringVector manifestNames_, arma::uvec personInSubset_){
  if(currentStatus < 3){
    Rcpp::stop("Please define the model matrices and add parameters as well as add derivative elements before calling addRawData.");
  }
  rawData = rawData_;
  personInSubset = personInSubset_;
  manifestNames = manifestNames_;
  currentStatus = 4;
  wasFit = false; // reset fit 
  return;
}

void SEMCpp::addSubset(int N_,
                       int observed_, // number of variables without missings
                       arma::uvec notMissing_, // vector with indices of non-missing values
                       // the following elements are only relevant for N>1
                       arma::mat covariance_,
                       arma::colvec means_,
                       // raw data is required for N == 1
                       arma::colvec rawData_){
  
  if(currentStatus < 4){
    Rcpp::stop("Please define the model matrices, add parameters as well as derivative elements and raw data before calling addSubset.");
  }
  
  data.addSubset(N_,
                 observed_, // number of variables without missings
                 notMissing_, // vector with indices of non-missing values
                 // the following elements are only relevant for N>1
                 covariance_,
                 means_,
                 // raw data is required for N == 1
                 rawData_);
  currentStatus = 5;
  wasFit = false; // reset fit 
  return;
  
}

void SEMCpp::removeSubset(int whichSubset){
  data.removeSubset(whichSubset);
  wasFit = false; // reset fit 
  return;
}

void SEMCpp::setMatrix(std::string whichMatrix, arma::mat values){
  currentStatus = 1;
  wasFit = false; // reset fit 
  
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
  wasFit = false; // reset fit 
  
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
                                  arma::vec rawValue_){
  if(currentStatus < 1){
    Rcpp::stop("Please define the model matrices before adding parameters");
  }
  parameterTable.initialize(label_,
                            location_,
                            row_,
                            col_,
                            value_,
                            rawValue_);
  currentStatus = 2;
  wasFit = false; // reset fit 
  return;
}

void SEMCpp::setParameters(Rcpp::StringVector label_,
                           arma::vec value_,
                           bool raw){
  wasFit = false; // reset fit 
  // step one: change parameters in parameterTable
  parameterTable.setParameters(label_, value_, raw);
  // step two: change parameters in model matrices
  std::string parLocation;
  for(int par = 0; par < parameterTable.label.length(); par++){
    parLocation = Rcpp::as<std::string>(parameterTable.location.at(par));
    if(parLocation.compare("Amatrix") == 0){
      Amatrix(parameterTable.row.at(par), 
              parameterTable.col.at(par)) = parameterTable.value.at(par);
      continue;
    }
    if(parLocation.compare("Smatrix") == 0){
      Smatrix(parameterTable.row.at(par), 
              parameterTable.col.at(par)) = parameterTable.value.at(par);
      continue;
    }
    if(parLocation.compare("Mvector") == 0){
      Mvector(parameterTable.row.at(par), 
              parameterTable.col.at(par)) = parameterTable.value.at(par);
      continue;
    }
    Rcpp::Rcout << parameterTable.label.at(par) << std::endl;
    Rcpp::stop("Could not find parameter in SEMCpp::setParameters");
  }
  
  // reset all fit elements
  m2LL = 0.0;
  impliedCovariance.fill(arma::datum::nan);
  impliedMeans.fill(arma::datum::nan);
  return;
}


void SEMCpp::addDerivativeElement(std::string label_, 
                                  std::string location_, 
                                  bool isVariance_, 
                                  arma::mat positionMatrix_){
  if(currentStatus < 2){
    Rcpp::stop("Please define the model matrices and add parameters before calling addDerivativeElement.");
  }
  derivElements.addDerivativeElement(label_, 
                                     location_, 
                                     isVariance_, 
                                     positionMatrix_);
  currentStatus = 3;
  wasFit = false; // reset fit 
  return;
}

// getter
Rcpp::DataFrame SEMCpp::getParameters(){
  
  return(parameterTable.getParameters());
  
}

Rcpp::StringVector SEMCpp::getParameterLabels(){
  Rcpp::StringVector uniqueParLabels(derivElements.uniqueLabels.size());
  for(int i = 0; i < derivElements.uniqueLabels.size(); i++){
    uniqueParLabels(i) = derivElements.uniqueLabels.at(i);
  }
  return(uniqueParLabels);
  
}
// fit functions


// [[Rcpp :: depends ( RcppArmadillo )]]

void SEMCpp::implied(){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  // if(!wasFit){
  //   Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  // }
  // step one: compute implied mean and covariance
  impliedCovariance = computeImpliedCovariance(Fmatrix, Amatrix, Smatrix);
  impliedMeans = computeImpliedMeans(Fmatrix, Amatrix, Mvector);
  
  return;
}


double SEMCpp::fit(){
  if(!wasChecked){
    wasChecked = checkModel();
  }
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
  if(!wasFit){
    Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  }
  arma::mat scoresMat = scores(*this, raw);
  
  return(scoresMat);
}

arma::rowvec SEMCpp::getGradients(bool raw){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  if(!wasFit){
    Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  }
  arma::rowvec gradients = gradientsByGroup(*this, raw);
  
  return(gradients);
}

arma::mat SEMCpp::getHessian(Rcpp::StringVector label_,
                             arma::vec value_,
                             bool raw,
                             double eps){
  if(!wasChecked){
    wasChecked = checkModel();
  }
  if(!wasFit){
    Rcpp::stop("The model has not been fitted yet. Call Model$fit() first.");
  }
  arma::mat hessian = approximateHessian(*this, 
                                         label_,
                                         value_,
                                         raw,
                                         eps);
  
  return(hessian);
}

RCPP_EXPOSED_CLASS(SEMCpp)
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
  }

