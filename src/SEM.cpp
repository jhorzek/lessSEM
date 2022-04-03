#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "implied.hpp"
#include "fit.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

void SEMCpp::addRawData(arma::mat rawData_, Rcpp::CharacterVector manifestNames_){
  rawData = rawData_;
  manifestNames = manifestNames_;
}

void SEMCpp::addSubset(int N_,
                    int observed_, // number of variables without missings
                    arma::uvec notMissing_, // vector with indices of non-missing values
                    // the following elements are only relevant for N>1
                    arma::mat covariance_,
                    arma::colvec means_,
                    // raw data is required for N == 1
                    arma::colvec dataNoMissing_){
  data.addSubset(N_,
                 observed_, // number of variables without missings
                 notMissing_, // vector with indices of non-missing values
                 // the following elements are only relevant for N>1
                 covariance_,
                 means_,
                 // raw data is required for N == 1
                 dataNoMissing_);
  return;
  
}

void SEMCpp::removeSubset(int whichSubset){
  data.removeSubset(whichSubset);
  return;
}

void SEMCpp::setMatrix(std::string whichMatrix, arma::mat values){

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
    mVector = values;
    return;
  }
  
  Rcpp::stop("Unknown whichVector. Can only assign m.");
}

void SEMCpp::initializeParameters(Rcpp::CharacterVector label_,
                               Rcpp::CharacterVector location_,
                               arma::uvec row_,
                               arma::uvec col_,
                               arma::vec value_,
                               arma::vec rawValue_){
  parameterTable.initialize(label_,
                            location_,
                            row_,
                            col_,
                            value_,
                            rawValue_);
  return;
}

void SEMCpp::setParameters(Rcpp::CharacterVector label_,
                        arma::vec value_,
                        bool raw){
  // step one: change parameters in parameterTable
  parameterTable.setParameters(label_, value_, raw);
  // step two: change parameters in model matrices
  
  for(int par = 0; par < parameterTable.label.length(); par++){
    
    if(parameterTable.label.at(par) == "Amatrix"){
      Amatrix(parameterTable.row.at(par), 
              parameterTable.col.at(par)) = parameterTable.value.at(par);
      continue;
    }
    if(parameterTable.label.at(par) == "Smatrix"){
      Smatrix(parameterTable.row.at(par), 
              parameterTable.col.at(par)) = parameterTable.value.at(par);
      continue;
    }
    if(parameterTable.label.at(par) == "mVector"){
      mVector(parameterTable.row.at(par), 
              parameterTable.col.at(par)) = parameterTable.value.at(par);
      continue;
    }
    Rcpp::stop("Could not find parameter in SEMCpp::setParameters");
  }
  
  // reset all fit elements
  m2LL = 0.0;
  impliedCovariance.fill(arma::datum::nan);
  impliedMeans.fill(arma::datum::nan);
  return;
}


// getter
Rcpp::DataFrame SEMCpp::getParameters(){
  
  return(parameterTable.getParameters());
  
}

// fit functions


// [[Rcpp :: depends ( RcppArmadillo )]]

void SEMCpp::implied(){
  
  // step one: compute implied mean and covariance
  impliedCovariance = computeImpliedCovariance(Fmatrix, Amatrix, Smatrix);
  impliedMeans = computeImpliedMeans(Fmatrix, Amatrix, mVector);
  
  return;
}


double SEMCpp::fit(){
  m2LL = 0.0;

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
                                                 currentSubset.dataNoMissing, 
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

RCPP_EXPOSED_CLASS(SEMCpp)
  RCPP_MODULE(SEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<SEMCpp>( "SEMCpp" )
      .constructor("Creates a new SEMCpp.")
      .field_readonly( "A", &SEMCpp::Amatrix, "Matrix with directed effects")
      .field_readonly( "S", &SEMCpp::Smatrix, "Matrix with undirected paths")
      .field_readonly( "F", &SEMCpp::Fmatrix, "Filter matrix to separate latent and observed variables")
      .field_readonly( "m", &SEMCpp::mVector, "Vector with means of observed and latent variables")
      .field_readonly( "impliedCovariance", &SEMCpp::impliedCovariance, "implied covariance matrix")
      .field_readonly( "impliedMeans", &SEMCpp::impliedMeans, "implied means vector")
      .field_readonly( "m2LL", &SEMCpp::m2LL, "minus 2 log-likelihood")
      .field_readonly( "rawData", &SEMCpp::rawData, "raw data set")
      .field_readonly( "manifestNames", &SEMCpp::manifestNames, "names of manifest variables")
    
    // methods
    .method( "setMatrix", &SEMCpp::setMatrix, "Fills the elements of a model matrix. Expects a char (A, S, or F), and a matrix with values")
    .method( "setVector", &SEMCpp::setVector, "Fills the elements of a model vector. Expects a char (m), and a vector with values")
    .method( "initializeParameters", &SEMCpp::initializeParameters, "Initializes the parameters of the model. Expects CharacterVector with labels, CharacterVector with location, uvec with rows, uvec with columns, vec with values, and vec with rawValues")
    .method( "setParameters", &SEMCpp::setParameters, "Changes the parameters of a model. Expects CharacterVector with labels, vec with values, and boolean indicating if the values are raw (TRUE) or transformed (FALSE).")
    .method( "addRawData", &SEMCpp::addRawData, "Adds a raw data set. Expects matrix and vector with labels of manifest variables.")
    .method( "addSubset", &SEMCpp::addSubset, "Adds a subset to the data. Expects the sample size N of the subset, the number of observed variables without missings, uvec with indices of notMissing values, matrix with covariance, colvec with means, raw data without missings.")
    .method( "implied", &SEMCpp::implied, "Computes implied means and covariance matrix")
    .method( "fit", &SEMCpp::fit, "Fits the model. Returns -2 log likelihood")
    .method( "getParameters", &SEMCpp::getParameters, "Returns a data frame with model parameters.")
    ;
  }

