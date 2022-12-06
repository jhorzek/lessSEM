# include "mgSEM.h"

void mgParameters::addTransformation(
    std::vector<bool> isTransformation_,
    SEXP transformationFunctionSEXP,
    Rcpp::List transformationList_)
{
  hasTransformations = true;
  isTransformation = isTransformation_;
  // create pointer to the transformation function 
  Rcpp::XPtr<transformationFunctionPtr> xpTransformationFunction(transformationFunctionSEXP);
  transformationFunction = *xpTransformationFunction; 
  transformationList = transformationList_;
}

void mgParameters::transform()
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

arma::mat mgParameters::getTransformationGradients(){
  if(!hasTransformations) Rcpp::stop("Does not have transformations.");
  int nRealParameters; // number of parameters that are not transformations 
  // of other parameters
  for (int i = 0; i < isTransformation.size(); i++){
    nRealParameters += !isTransformation.at(i);
  }
  arma::mat currentGradients(uniqueLabels.size(), 
                             nRealParameters,
                             arma::fill::zeros);
  Rcpp::NumericVector parameterValues(uniqueLabels.size());
  Rcpp::CharacterVector parameterLabelsRcpp(uniqueLabels.size());
  arma::uvec selectRows(uniqueLabels.size());
  arma::colvec stepForward, stepBackward;
  std::string currentParameter;
  int j = 0;
  for(int i = 0; i < isTransformation.size(); i++){
    parameterValues.at(i) = uniqueValues.at(i);
    
    if(isTransformation.at(i)){
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
    if(isTransformation.at(i)) continue;
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


//' @name mgSEM
//' 
//' @title mgSEM class
//' 
//' @description internal mgSEM representation
//' 
//' @field new Creates a new mgSEM.
//' @field addModel add a model. Expects Rcpp::List
//' @field addTransformation adds transforamtions to a model
//' @field implied Computes implied means and covariance matrix
//' @field fit Fits the model. Returns -2 log likelihood
//' @field getParameters Returns a data frame with model parameters.
//' @field getParameterLabels Returns a vector with unique parameter labels as used internally.
//' @field getGradients Returns a matrix with scores.
//' @field getScores Returns a matrix with scores. Not yet implemented
//' @field getHessian Returns the hessian of the model. Expects the labels of the 
//' parameters and the values of the parameters as well as a boolean indicating if 
//' these are raw. Finally, a double (eps) controls the precision of the approximation.
//' @field computeTransformations compute the transformations.
//' 
void mgSEM::addModel(Rcpp::List SEMList){
  if(parameters.hasTransformations){
    Rcpp::stop("It seems like transformations were already added to the model. You cannot add further models.");
  }
  SEMCpp newModel;
  newModel.fill(SEMList);
  if(newModel.hasTransformations){
    Rcpp::stop("There should be not transformations in the sub-models.");
  }
  models.push_back(newModel);
  
  Rcpp::DataFrame newParameters = newModel.getParameters();
  Rcpp::StringVector newLabels = newParameters["label"];
  Rcpp::NumericVector newValues = newParameters["value"];
  
  parameters.modelParameterValues.push_back(newValues);
  parameters.modelParameterLabels.push_back(newLabels);
  
  parameters.modelParameterValues.at(models.size() - 1).names() = newLabels;
  
  std::string parameterLabel;
  
  // for each parameter, we check if this parameter is already present in
  // one of the existing models. If not, we add the 
  // parameter to our parameter vector.
  for(int i = 0; i < newLabels.length(); i++){
    
    parameterLabel = Rcpp::as< std::string >(newLabels.at(i));
    
    if(findStringInVector(parameterLabel, parameters.uniqueLabels, false) == -1){
      // not yet in the parameter vector
      parameters.uniqueValues.resize(parameters.uniqueValues.size()+1);
      parameters.uniqueValues.at(parameters.uniqueValues.size()-1) = newValues.at(i);
      
      parameters.uniqueLabelsRcpp.push_back(newLabels.at(i));
      parameters.uniqueLabels.push_back(Rcpp::as<std::string>(newLabels.at(i)));
      
    }
  } // end for
  
  // next, we update the indices for Rcpp and armadillo vectors
  parameters.parameterLocationInModelRcpp.resize(models.size());
  parameters.parameterLocationInVectorRcpp.resize(models.size());
  parameters.parameterLocationInModelUvec.resize(models.size());
  parameters.parameterLocationInVectorUvec.resize(models.size());
  
  for(int m = 0; m < models.size(); m++){
    // initialize as empty vector:
    Rcpp::IntegerVector currentInteger = {};
    parameters.parameterLocationInModelRcpp.at(m) = currentInteger;
    parameters.parameterLocationInVectorRcpp.at(m) = currentInteger;
    
    Rcpp::StringVector currentLabels = models.at(m).getParameterLabels();
    
    // now, let's fill the vectors:
    for(int p = 0; p < parameters.uniqueLabels.size(); p++){
      // check if this parameter is in the current model; if so, add its position 
      // to the indices
      int locatedInVector = findStringInVector(parameters.uniqueLabels.at(p), currentLabels, false);
      if(locatedInVector != -1){
        parameters.parameterLocationInModelRcpp.at(m).push_back(locatedInVector);
        parameters.parameterLocationInVectorRcpp.at(m).push_back(p);
      }
    }
    // also replace the uvecs for Rcpp:
    parameters.parameterLocationInModelUvec.at(m) = Rcpp::as<arma::uvec>(parameters.parameterLocationInModelRcpp.at(m));
    parameters.parameterLocationInVectorUvec.at(m) = Rcpp::as<arma::uvec>(parameters.parameterLocationInVectorRcpp.at(m));
  }
  
  // update length of gradients vector
  parameters.uniqueGradients.set_size(parameters.uniqueLabels.size());
  parameters.uniqueGradients.fill(arma::fill::zeros);
  
  parameters.uniqueHessian.set_size(parameters.uniqueLabels.size(), parameters.uniqueLabels.size());
  parameters.uniqueHessian.fill(arma::fill::zeros);
  
}// end addModel

void mgSEM::addTransformation(std::vector<bool> isTransformation_,
                              SEXP transformationFunctionSEXP,
                              Rcpp::List transformationList){
  parameters.addTransformation(isTransformation_, 
                               transformationFunctionSEXP,
                               transformationList
                               );
}

void mgSEM::computeTransformations(){
  parameters.transform();
}

// changing parameters

void mgSEM::setParameters(Rcpp::StringVector label_,
                          arma::vec value_,
                          bool raw){
  // change the global parameters
  int loc;
  for(int i = 0; i < label_.size(); i++){
    // change parameter
    loc = findStringInVector(Rcpp::as<std::string>(label_.at(i)), parameters.uniqueLabels, true);
    
    parameters.uniqueValues.at(loc) = value_.at(i);
    
  }
  
  // compute the transformations
  Rcpp::warning("Transformations missing");
  computeTransformations();
  
  // update the parameters in the models
  
  // from the map we can directly see which parameters should be changed in the 
  // models. We first change each model's parameter vector:
  for(int m = 0; m < models.size(); m++){
    models.at(m).setParameters(parameters.uniqueLabelsRcpp[parameters.parameterLocationInVectorRcpp.at(m)],
              parameters.uniqueValues.elem(parameters.parameterLocationInVectorUvec.at(m)),
              true
    );
  }
} // end set parameters

// getter
Rcpp::NumericVector mgSEM::getParameters(){
  Rcpp::NumericVector param(parameters.uniqueValues.size());
  for(int p = 0; p < param.size(); p++){
    param.at(p) = parameters.uniqueValues.at(p);
  }
  param.names() = parameters.uniqueLabelsRcpp;
  return(param);
}

Rcpp::StringVector mgSEM::getParameterLabels(){
  return(parameters.uniqueLabelsRcpp);
}

// fit related functions
void mgSEM::implied(){
  // compute implied means and covariance for each model
  for(int m = 0; m < models.size(); m++){
    models.at(m).implied();
  }
}

double mgSEM::fit(){
  
  double m2LL = 0.0;
  // compute fit for each model
  for(int m = 0; m < models.size(); m++){
    m2LL += models.at(m).fit();
  }
  
  return(m2LL);
}

arma::rowvec mgSEM::getGradients(bool t){
  arma::rowvec modelGradients;
  arma::mat transformationGradients;
  parameters.uniqueGradients.fill(arma::fill::zeros);
  
  // compute gradients for each model
  for(int m = 0; m < models.size(); m++){
    
    modelGradients = models.at(m).getGradients(true);
    
    // add the models gradients to the existing gradients:
    parameters.uniqueGradients.elem(parameters.parameterLocationInModelUvec.at(m)) += 
      modelGradients.elem(parameters.parameterLocationInModelUvec.at(m));
  }
  
  if(!parameters.hasTransformations){
    return(parameters.uniqueGradients);
  }
  
  // compute the transformation gradients
  transformationGradients = parameters.getTransformationGradients();
  
  return(parameters.uniqueGradients*transformationGradients);
}

arma::mat getScores(){
  Rcpp::stop("Score function not yet impemented for multi-group models");
  arma::mat m;
  return(m);
}

arma::mat mgSEM::getHessian(Rcpp::StringVector label_,
                            arma::vec value_,
                            double eps){
  
  Hessian = approximateHessian(*this, 
                               label_,
                               value_,
                               true,
                               eps);
  
  return(Hessian);
}


RCPP_MODULE(mgSEM_Cpp){
  using namespace Rcpp;
  Rcpp::class_<mgSEM>( "mgSEM" )
    .constructor("Creates a new SEMCpp.")
  // methods
  .method( "addModel", &mgSEM::addModel, "Adds a model. Expects and Rcpp::List")
  .method( "implied", &mgSEM::implied, "Computes implied means and covariance matrix")
  .method( "fit", &mgSEM::fit, "Fits the model. Returns -2 log likelihood")
  .method( "setParameters", &mgSEM::setParameters, "Set the parameters of a model.")
  .method( "getParameters", &mgSEM::getParameters, "Returns a data frame with model parameters.")
  .method( "getParameterLabels", &mgSEM::getParameterLabels, "Returns a vector with unique parameter labels as used internally.")
  .method( "getGradients", &mgSEM::getGradients, "Returns a matrix with scores.")
  //.method( "getScores", &mgSEM::getScores, "Returns a matrix with scores.")
  .method( "getHessian", &mgSEM::getHessian, "Returns the hessian of the model. Expects the labels of the parameters and the values of the parameters as well as a boolean indicating if these are raw. Finally, a double (eps) controls the precision of the approximation.")
  .method( "addTransformation", &mgSEM::addTransformation, "Add a transformation function. Expects parameterLabels and pointer to function.")
  .method( "computeTransformations", &mgSEM::computeTransformations, "Compute all transformations")
  ;
}