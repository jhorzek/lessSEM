#ifndef MGSEMMODULE_H
#define MGSEMMODULE_H

#include <RcppArmadillo.h>
#include "SEM.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

int findStringInVector(std::string what, std::vector<std::string> where, bool throwError){
  
  for (int i = 0; i < where.size(); i++){
    if(where.at(i).compare(what) == 0){
      return(i);
    }
  }
  if(throwError){
    Rcpp::stop("Could not find parameter.");
  }else{
    return -1;
  }
}

class mgParameters{
public:
  
  Rcpp::NumericVector uniqueValues;
  std::vector<std::string> uniqueLabels;
  arma::colvec uniqueGradients;
  arma::mat uniqueHessian;
  std::vector<bool> isTransformation;
  bool hasTransformations = false;
  
  // IDEA: REPLACE THE FOLLOWING WITH A BOOLEAN MATRIX INDICATING FOR
  // EACH PARAMETER IF THIS PARAMETER IS IN A SPECIFIC MODEL. THIS
  // ALLOWS FOR SIMPLE SUBSETTING OF THE VECTOR WHICH IS USEFUL WHEN COMPUTING
  // THE GRADIENTS
  std::map< std::string, std::vector<int>> parameterInModel; // vector telling us in which model a parameter is located
  
  std::vector<std::vector<int>> model; // saves the models where the parameters are located
  
  // vectors with values and labels for each model
  std::vector<Rcpp::NumericVector> modelParameterValues;
  std::vector<Rcpp::StringVector> modelParameterLabels;
  
  // in case of transformations
  transformationFunctionPtr transformationFunction;
  Rcpp::List transformationList;
  
  // we use the same initialize function as that of parameters.
  // Here, we initialize all raw and transformed parameters
  
  // Now, we have to tell our SEM, which parameters are transformations
  // and how to compute those
  void addTransformation(SEXP transformationFunctionSEXP,
                         Rcpp::List transformationList_)
    {
    hasTransformations = true;
    // create pointer to the transformation function 
    Rcpp::XPtr<transformationFunctionPtr> xpTransformationFunction(transformationFunctionSEXP);
    transformationFunction = *xpTransformationFunction; 
    transformationList = transformationList_;
  }
  
  void transform()
  {
    Rcpp::NumericVector params(uniqueLabels.size());
    Rcpp::CharacterVector paramLabels(uniqueLabels.size());
    for(int i = 0; i < uniqueLabels.size(); i++){
      params.at(i) = uniqueValues.at(i);
      paramLabels.at(i) = uniqueLabels.at(i);
    }
    params.names() = paramLabels;
    
    params = transformationFunction(params, transformationList);
    
    // also change the parameter values in the unique parameter vector. 
    // these are the ones that are actually used internally
    std::string parameterLabel;
    int location;
    for(int p = 0; p < paramLabels.length(); p++){
      parameterLabel = Rcpp::as< std::string >(paramLabels.at(p));
      location = findStringInVector(parameterLabel, uniqueLabels, true);
      uniqueValues.at(location) = params.at(p);
    }
  }
  
  arma::mat getTransformationGradients(){
    if(!hasTransformations) Rcpp::stop("Does not have transformations.");
    arma::mat currentGradients(nModelParameters, 
                               nTransformationParameters,
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
    
    double eps = 1e-6;
    
    j = 0;
    for(int i = 0; i < parameterValues.size(); i++){
      currentParameter = parameterLabelsRcpp.at(i);
      if(parameterMap.at(currentParameter).isTransformation) continue;
      parameterValues.at(i) += eps;
      
      stepForward = Rcpp::as<arma::colvec>(transformationFunction(parameterValues, transformationList)).rows(selectRows);
      
      parameterValues.at(i) -= 2.0*eps;
      
      stepBackward = Rcpp::as<arma::colvec>(transformationFunction(parameterValues, transformationList)).rows(selectRows);
      
      parameterValues.at(i) += eps;
      
      currentGradients.col(j) = (stepForward - stepBackward)/(2.0*eps);
      j++;
    }
    
    return(currentGradients);
  }
  
};



class mgSEM{
public:
  std::vector<SEMCpp> models;
  
  mgParameters parameters;
  
  Rcpp::StringVector parameterLabels;
  Rcpp::NumericVector parameterValues;
  std::vector<int> locatedInModel;
  
  void addModel(SEMCpp newModel){
    models.push_back(newModel);
    
    Rcpp::DataFrame newParameters = newModel.getParameters();
    Rcpp::StringVector newLabels = newParameters["label"];
    Rcpp::NumericVector newValues = newParameters["value"];
    
    parameters.modelParameterValues.push_back(newValues);
    parameters.modelParameterLabels.push_back(newLabels);
    
    parameters.modelParameterValues.at(models.size() - 1).names() = newLabels;
    
    std::string parameterLabel;
    
    for(int i = 0; i < newLabels.length(); i++){
      
      parameterLabel = Rcpp::as< std::string >(newLabels.at(i));
      
      if(findStringInVector(parameterLabel, parameters.uniqueLabels, false) == -1){
        // not yet in the parameter vector
        parameters.uniqueValues.push_back(newValues.at(i));
        parameters.uniqueLabels.push_back(Rcpp::as<std::string>(newLabels.at(i)));
        parameters.parameterInModel[parameterLabel].push_back(models.size()-1); // indices starting with 0
        
      }else{
        
        parameters.parameterInModel[parameterLabel].push_back(models.size()-1); // indices starting with 0
        
      }
      
    } // end for
    
    // update names of parameters
    parameters.uniqueValues.names() = parameters.uniqueLabels;
    // update length of gradients vector
    parameters.uniqueGradients.set_size(parameters.uniqueLabels.size());
    parameters.uniqueGradients.fill(arma::fill::zeros);
    
    parameters.uniqueHessian.set_size(parameters.uniqueLabels.size(), parameters.uniqueLabels.size());
    parameters.uniqueHessian.fill(arma::fill::zeros);
    
  }// end addModel
  
  void addTransformation(SEXP transformationFunctionSEXP,
                         Rcpp::List transformationList);
  
  void computeTransformations();
  
  // changing parameters
  
  void setParameters(Rcpp::StringVector label_,
                     arma::vec value_,
                     bool raw){
    // change the global parameters
    for(int i = 0; i < label_.size(); i++){
      // change parameter
      parameters.uniqueValues[Rcpp::as<std::string>(label_.at(i))] = value_.at(i);
      
    }
    
    // compute the transformations
    Rcpp::warning("Transformations missing");
    computeTransformations();
    
    // update the parameters in the models
    
    // from the map we can directly see which parameters should be changed in the 
    // models. We first change each model's parameter vector:
    for(int p = 0; p < parameters.uniqueLabels.size(); p++){
      for(int m = 0; m < parameters.parameterInModel[parameters.uniqueLabels.at(p)].size(); m++){
        parameters.modelParameterValues.at(m)[parameters.uniqueLabels.at(p)] = parameters.uniqueValues.at(p);
      }
    }
    
    // Having model-specific parameter vectors, we can now change the parameters of each model
    for(int m = 0; m < models.size(); m++){
      models.at(m).setParameters(parameters.modelParameterLabels.at(m),  
                parameters.modelParameterValues.at(m), 
                false);
    }
    // end set parameters
  }
  
  // getter
  Rcpp::DataFrame getParameters();
  Rcpp::StringVector getParameterLabels(){
    Rcpp::StringVector labels(parameters.uniqueLabels.size());
    for(int i = 0; i < parameters.uniqueLabels.size(); i++){
      labels.at(i) = parameters.uniqueLabels.at(i);
    }
    return(labels);
  }
  
  // fit related functions
  void implied(){
    // compute implied means and covariance for each model
    for(int m = 0; m < models.size(); m++){
      models.at(m).implied();
    }
  }
  
  double fit(){
    
    double m2LL = 0.0;
    // compute fit for each model
    for(int m = 0; m < models.size(); m++){
      m2LL += models.at(m).fit();
    }
    
    return(m2LL);
  }
  
  arma::rowvec getGradients(){
    
    arma::rowvec modelGradients;
    Rcpp::StringVector modelNames;
    parameters.uniqueGradients.fill(arma::fill::zeros);
    
    // compute gradients for each model
    for(int m = 0; m < models.size(); m++){
      
      modelGradients = models.at(m).getGradients(true);
      modelNames = models.at(m).getParameterLabels();
      
      // add the models gradients to the existing gradients:
      for(int p = 0; p < modelNames.size(); p ++){
        if(parameters.uniqueLabels.at(p).compare(Rcpp::as<std::string>(modelNames.at(p))) == 0){
          parameters.uniqueGradients.col(p) += modelGradients.col(p);
        }
      }
    }
    
    return(parameters.uniqueGradients);
  }
  
  arma::mat getScores(){
    Rcpp::stop("Score function not yet impemented for multi-group models");
  }
  
  arma::mat getHessian(Rcpp::StringVector label_,
                       arma::vec value_,
                       double eps){
    setParameters(label_, value_, true);

    std::string par1, par2;
    int rowAt, colAt;
    arma::mat modelHessian;
    Rcpp::StringVector modelLabels;
    
    for(int m = 0; m < models.size(); m++){
      
      modelHessian = models.at(m).getHessian(parameters.modelParameterLabels.at(m),
                               parameters.modelParameterValues.at(m),
                               false,
                               eps
                               );
      modelLabels = models.at(m).getParameterLabels();
      
      for(int r = 0; r < modelHessian.n_rows; r++){
        par1 = Rcpp::as<std::string>(modelLabels.at(r));
        rowAt = findStringInVector(par1, parameters.uniqueLabels, true);
        
        for(int c = 0; c < modelHessian.n_cols; c++){
          par2 = Rcpp::as<std::string>(modelLabels.at(c));
          colAt = findStringInVector(par2, parameters.uniqueLabels, true);
          parameters.uniqueHessian(rowAt, colAt) += modelHessian(r,c);
        }
      }
    }// end for model
    
    return(parameters.uniqueHessian);
  }
};

# endif