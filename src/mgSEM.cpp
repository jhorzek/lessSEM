# include "mgSEM.h"

void mgParameters::addTransformation(
    Rcpp::NumericVector extendedParameters,
    std::vector<bool> isTransformation_,
    SEXP transformationFunctionSEXP,
    Rcpp::List transformationList_)
{
  // first, let's check that all parameters have been added to the END of the
  // existing parameter vector
  Rcpp::StringVector extendedLabels = extendedParameters.names();
  
  for(unsigned int p = 0; p < uniqueLabelsRcpp.size(); p++){
    
    if(uniqueLabelsRcpp.at(p) != extendedLabels.at(p)) Rcpp::stop("Mismatch in parameters");
    
  }
  
  // extend existing vectors
  uniqueLabelsRcpp = extendedLabels;
  for(unsigned int p = uniqueLabels.size(); p < extendedLabels.size(); p++){
    uniqueLabels.push_back(Rcpp::as<std::string>(extendedLabels.at(p)));
  }
  uniqueGradients.resize(extendedLabels.size());
  uniqueHessian.resize(extendedLabels.size(), extendedLabels.size());
  uniqueValues.resize(extendedLabels.size());
  for(unsigned int p = 0; p < extendedParameters.size(); p++){
    uniqueValues.at(p) = extendedParameters.at(p);
  }
  
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
  for(unsigned int i = 0; i < uniqueLabels.size(); i++){
    params.at(i) = uniqueValues.at(i);
    paramLabels.at(i) = uniqueLabels.at(i);
  }
  params.names() = paramLabels;
  
  params = transformationFunction(params, transformationList);
  
  // also change the parameter values in the unique parameter vector. 
  // these are the ones that are actually used internally
  std::string parameterLabel;
  int location;
  for(unsigned int p = 0; p < paramLabels.length(); p++){
    parameterLabel = Rcpp::as< std::string >(paramLabels.at(p));
    location = findStringInVector(parameterLabel, uniqueLabels, true);
    uniqueValues.at(location) = params.at(p);
  }
}

arma::mat mgParameters::getTransformationGradients(){
  if(!hasTransformations) Rcpp::stop("Does not have transformations.");
  int nRealParameters = 0; // number of parameters that are not transformations 
  // of other parameters
  for(unsigned int i = 0; i < isTransformation.size(); i++){
    nRealParameters += !isTransformation.at(i);
  }
  arma::mat currentGradients(uniqueLabels.size(), 
                             nRealParameters,
                             arma::fill::zeros);
  Rcpp::NumericVector parameterValues(uniqueLabels.size());
  arma::uvec selectRows(nRealParameters);
  arma::colvec stepForward, stepBackward;
  int j = 0;
  for(unsigned int i = 0; i < isTransformation.size(); i++){
    parameterValues.at(i) = uniqueValues.at(i);
    if(isTransformation.at(i)){
      // we want to remove all parameters which are only in transformations and not in the model
      continue;
    }
    //selectRows(j) = i;
    j++;
  }
  parameterValues.names() = uniqueLabels;
  
  j = 0;
  for(unsigned int i = 0; i < parameterValues.size(); i++){
    if(isTransformation.at(i)) continue;
    parameterValues.at(i) += gradientStepSize;
    
    stepForward = Rcpp::as<arma::colvec>(transformationFunction(parameterValues, transformationList));//.rows(selectRows);
    
    parameterValues.at(i) -= 2.0*gradientStepSize;
    
    stepBackward = Rcpp::as<arma::colvec>(transformationFunction(parameterValues, transformationList));//.rows(selectRows);
    
    parameterValues.at(i) += gradientStepSize;
    
    currentGradients.col(j) = (stepForward - stepBackward)/(2.0*gradientStepSize);
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
//' @field setTransformationGradientStepSize change the step size of the gradient computation for the transformations
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
  sampleSize += newModel.sampleSize;
  
  Rcpp::DataFrame newParameters = newModel.getParameters();
  Rcpp::StringVector newLabels = newParameters["label"];
  Rcpp::NumericVector newValues = newParameters["rawValue"];
  
  std::string parameterLabel;
  
  // for each parameter, we check if this parameter is already present in
  // one of the existing models. If not, we add the 
  // parameter to our parameter vector.
  for(unsigned int i = 0; i < newLabels.length(); i++){
    
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
  
  for(unsigned int m = 0; m < models.size(); m++){
    // initialize as empty vector:
    Rcpp::IntegerVector currentInteger = {};
    parameters.parameterLocationInModelRcpp.at(m) = currentInteger;
    parameters.parameterLocationInVectorRcpp.at(m) = currentInteger;
    
    Rcpp::StringVector currentLabels = models.at(m).getParameterLabels();
    
    // now, let's fill the vectors:
    for(unsigned int p = 0; p < parameters.uniqueLabels.size(); p++){
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

void mgSEM::addTransformation(Rcpp::NumericVector extendedParameters,
                              std::vector<bool> isTransformation_,
                              SEXP transformationFunctionSEXP,
                              Rcpp::List transformationList_){
  parameters.addTransformation(extendedParameters,
                               isTransformation_,
                               transformationFunctionSEXP,
                               transformationList_);
}

void mgSEM::computeTransformations(){
  parameters.transform();
}

// changing parameters

void mgSEM::setParameters(Rcpp::StringVector label_,
                          arma::vec value_,
                          bool raw){
  if(!raw) Rcpp::stop("Cannot set parameters for non-raw values");
  // change the global parameters
  int loc;
  for(unsigned int i = 0; i < label_.size(); i++){
    // change parameter
    loc = findStringInVector(Rcpp::as<std::string>(label_.at(i)), parameters.uniqueLabels, true);
    
    parameters.uniqueValues.at(loc) = value_.at(i);
    
  }
  
  // compute the transformations
  if(parameters.hasTransformations) computeTransformations();
  
  // update the parameters in the models
  for(unsigned int m = 0; m < models.size(); m++){
    models.at(m).setParameters(parameters.uniqueLabelsRcpp[parameters.parameterLocationInVectorRcpp.at(m)],
              parameters.uniqueValues.elem(parameters.parameterLocationInVectorUvec.at(m)),
              true
    );
  }
} // end set parameters

// getter
Rcpp::List mgSEM::getParameters(){
  Rcpp::NumericVector param(parameters.uniqueValues.size());
  for(unsigned int p = 0; p < param.size(); p++){
    param.at(p) = parameters.uniqueValues.at(p);
  }
  param.names() = parameters.uniqueLabelsRcpp;
  Rcpp::List retPar = Rcpp::List::create(Rcpp::Named("parmeters") = param,
                                         Rcpp::Named("isTransformation") = parameters.isTransformation);
  return(retPar);
}

Rcpp::List mgSEM::getSubmodelParameters(){
  Rcpp::List paramList;
  for(unsigned int m = 0; m < models.size(); m++){
    paramList.push_back(models.at(m).getParameters());
  }
  return(paramList);
}

Rcpp::StringVector mgSEM::getParameterLabels(){
  return(parameters.uniqueLabelsRcpp);
}

// fit related functions
void mgSEM::implied(){
  // compute implied means and covariance for each model
  for(unsigned int m = 0; m < models.size(); m++){
    models.at(m).implied();
  }
}

bool mgSEM::impliedIsPD(){
  bool isPd = true;
  for(unsigned int m = 0; m < models.size(); m++){
    isPd = isPd && models.at(m).impliedIsPD();
  }
  return(isPd);
}

double mgSEM::fit(){
  
  m2LL = 0.0;
  // compute fit for each model
  for(unsigned int m = 0; m < models.size(); m++){
    m2LL += models.at(m).fit();
  }
  
  return(m2LL);
}

arma::rowvec mgSEM::getGradients(bool raw){
  if(!raw) Rcpp::stop("Cannot compute gradients for non-raw values.");
  arma::rowvec modelGradients;
  arma::mat transformationGradients;
  parameters.uniqueGradients.fill(arma::fill::zeros);
  
  // compute gradients for each model
  for(unsigned int m = 0; m < models.size(); m++){
    modelGradients = models.at(m).getGradients(true);
    // add the models gradients to the existing gradients:
    parameters.uniqueGradients.elem(parameters.parameterLocationInVectorUvec.at(m)) += 
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
                            bool raw,
                            double eps){
  if(!raw) Rcpp::stop("Cannot compute Hessian for non-raw values.");
  Hessian = approximateHessian(*this, 
                               label_,
                               value_,
                               true,
                               eps);
  
  return(Hessian);
}

void mgSEM::setTransformationGradientStepSize(double gradientStepSize){
  parameters.gradientStepSize = gradientStepSize;
}

RCPP_MODULE(mgSEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<mgSEM>( "mgSEM" )
    .constructor("Creates a new SEMCpp.")
  // fields
  .field_readonly("sampleSize", &mgSEM::sampleSize, "Sum of all N")
  .field_readonly("m2LL", &mgSEM::m2LL, "-2 log-Likelihood")
  // methods
  .method( "addModel", &mgSEM::addModel, "Adds a model. Expects and Rcpp::List")
  .method( "implied", &mgSEM::implied, "Computes implied means and covariance matrix")
  .method( "fit", &mgSEM::fit, "Fits the model. Returns -2 log likelihood")
  .method( "setParameters", &mgSEM::setParameters, "Set the parameters of a model.")
  .method( "getParameters", &mgSEM::getParameters, "Returns a vector with raw model parameters.")
  .method( "getSubmodelParameters", &mgSEM::getSubmodelParameters, "Returns a list with parameters for each model.")
  .method( "getParameterLabels", &mgSEM::getParameterLabels, "Returns a vector with unique parameter labels as used internally.")
  .method( "getGradients", &mgSEM::getGradients, "Returns a matrix with scores.")
  //.method( "getScores", &mgSEM::getScores, "Returns a matrix with scores.")
    .method( "getHessian", &mgSEM::getHessian, "Returns the hessian of the model. Expects the labels of the parameters and the values of the parameters as well as a boolean indicating if these are raw. Finally, a double (eps) controls the precision of the approximation.")
    .method( "addTransformation", &mgSEM::addTransformation, "Add a transformation function. Expects parameterLabels and pointer to function.")
    .method( "computeTransformations", &mgSEM::computeTransformations, "Compute all transformations")
    .method( "setTransformationGradientStepSize", &mgSEM::setTransformationGradientStepSize, "Change the step size used in the computation of the transformation gradients.")
  ;
}