#include <RcppArmadillo.h>
#include "SEM.h"
#include "gpFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

//'@name istaLspGeneralPurpose
//'@title lsp optimization with ista
//'@description Object for lsp optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'an R function to compute the fit, an R function to compute the gradients, a
//'list with elements the fit and gradient function require, a theta and a lambda value (alpha must be 1).
//'@returns a list with fit results
class istaLspGeneralPurpose{
  public:
    
    
    // settings
  const Rcpp::NumericVector weights;
  
  // control optimizer
  const double L0;
  const double eta;
  const bool accelerate;
  const int maxIterOut;
  const int maxIterIn;
  const double breakOuter;
  const lessSEM::convCritInnerIsta convCritInner;
  const double sigma;
  const lessSEM::stepSizeInheritance stepSizeInh;
  const int verbose; 
  
  // constructor
  istaLspGeneralPurpose(
    const Rcpp::NumericVector weights_,
    Rcpp::List control
  ):
    weights(weights_),
  L0(Rcpp::as<double> (control["L0"])),
  eta(Rcpp::as<double> (control["eta"])),
  accelerate(Rcpp::as<bool> (control["accelerate"])),
  maxIterOut(Rcpp::as<int> (control["maxIterOut"])),
  maxIterIn(Rcpp::as<int> (control["maxIterIn"])),
  breakOuter(Rcpp::as<double> (control["breakOuter"])),
  convCritInner(static_cast<lessSEM::convCritInnerIsta>(Rcpp::as<int> (control["convCritInner"]))),
  sigma(Rcpp::as<double> (control["sigma"])),
  stepSizeInh(static_cast<lessSEM::stepSizeInheritance>(Rcpp::as<int> (control["stepSizeInheritance"]))),
  verbose(Rcpp::as<int> (control["verbose"])){}
  
  Rcpp::List optimize(
    Rcpp::NumericVector startingValues_, 
    Rcpp::Function fitFunction,
    Rcpp::Function gradientFunction,
    Rcpp::List userSuppliedElements,
    double theta_,
    double lambda_){
    
    generalPurposeFitFramework gpFF(fitFunction, gradientFunction, userSuppliedElements);
    
    int sampleSize = 1;
    
    lessSEM::tuningParametersLSP tp;
    tp.theta = theta_;
    tp.lambda = lambda_;
    tp.weights = weights;
    
    // we won't need the smooth penalty; but we need to specify some tuning 
    // parameters for the function call
    lessSEM::tuningParametersLSP smoothTp;
    smoothTp.lambda = 0.0;
    
    lessSEM::proximalOperatorLSP proximalOperatorLSP_;
    lessSEM::penaltyLSP penalty_;
    lessSEM::noSmoothPenalty<lessSEM::tuningParametersLSP> smoothPenalty_;
    
    lessSEM::control controlIsta = {
      L0,
      eta,
      accelerate,
      maxIterOut,
      maxIterIn,
      breakOuter,
      convCritInner,
      sigma,
      stepSizeInh,
      sampleSize,
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::ista(
      gpFF,
      startingValues_,
      proximalOperatorLSP_,
      penalty_,
      smoothPenalty_,
      tp,
      smoothTp,
      controlIsta
    );
    
    Rcpp::NumericVector finalParameters(fitResults_.parameterValues.n_elem);
    for(int i = 0; i < fitResults_.parameterValues.n_elem; i++){
      finalParameters.at(i) = fitResults_.parameterValues.at(i);
    }
    finalParameters.names() = startingValues_.names();
    
    if(!fitResults_.convergence) Rcpp::warning("Optimizer did not converge");
    Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("fit") = fitResults_.fit,
      Rcpp::Named("convergence") = fitResults_.convergence,
      Rcpp::Named("rawParameters") = finalParameters,
      Rcpp::Named("fits") = fitResults_.fits
    );
    return(result);
  }
};

//'@name istaLspGeneralPurposeCpp
//'@title lsp optimization with ista
//'@description Object for lsp optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEXP function pointer to compute the fit, a SEXP function pointer to compute the gradients, a
//'list with elements the fit and gradient function require, a theta and a lambda value (alpha must be 1).
//'@returns a list with fit results
class istaLspGeneralPurposeCpp{
public:
  
  
  // settings
  const Rcpp::NumericVector weights;
  
  // control optimizer
  const double L0;
  const double eta;
  const bool accelerate;
  const int maxIterOut;
  const int maxIterIn;
  const double breakOuter;
  const lessSEM::convCritInnerIsta convCritInner;
  const double sigma;
  const lessSEM::stepSizeInheritance stepSizeInh;
  const int verbose; 
  
  // constructor
  istaLspGeneralPurposeCpp(
    const Rcpp::NumericVector weights_,
    Rcpp::List control
  ):
    weights(weights_),
    L0(Rcpp::as<double> (control["L0"])),
    eta(Rcpp::as<double> (control["eta"])),
    accelerate(Rcpp::as<bool> (control["accelerate"])),
    maxIterOut(Rcpp::as<int> (control["maxIterOut"])),
    maxIterIn(Rcpp::as<int> (control["maxIterIn"])),
    breakOuter(Rcpp::as<double> (control["breakOuter"])),
    convCritInner(static_cast<lessSEM::convCritInnerIsta>(Rcpp::as<int> (control["convCritInner"]))),
    sigma(Rcpp::as<double> (control["sigma"])),
    stepSizeInh(static_cast<lessSEM::stepSizeInheritance>(Rcpp::as<int> (control["stepSizeInheritance"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      SEXP fitFunction,
      SEXP gradientFunction,
      Rcpp::List userSuppliedElements,
      double theta_,
      double lambda_){
    
    generalPurposeFitFrameworkCpp gpFF(startingValues_, fitFunction, gradientFunction, userSuppliedElements);
    
    int sampleSize = 1;
    
    lessSEM::tuningParametersLSP tp;
    tp.theta = theta_;
    tp.lambda = lambda_;
    tp.weights = weights;
    
    // we won't need the smooth penalty; but we need to specify some tuning 
    // parameters for the function call
    lessSEM::tuningParametersLSP smoothTp;
    smoothTp.lambda = 0.0;
    
    lessSEM::proximalOperatorLSP proximalOperatorLSP_;
    lessSEM::penaltyLSP penalty_;
    lessSEM::noSmoothPenalty<lessSEM::tuningParametersLSP> smoothPenalty_;
    
    lessSEM::control controlIsta = {
      L0,
      eta,
      accelerate,
      maxIterOut,
      maxIterIn,
      breakOuter,
      convCritInner,
      sigma,
      stepSizeInh,
      sampleSize,
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::ista(
      gpFF,
      startingValues_,
      proximalOperatorLSP_,
      penalty_,
      smoothPenalty_,
      tp,
      smoothTp,
      controlIsta
    );
    
    Rcpp::NumericVector finalParameters(fitResults_.parameterValues.n_elem);
    for(int i = 0; i < fitResults_.parameterValues.n_elem; i++){
      finalParameters.at(i) = fitResults_.parameterValues.at(i);
    }
    finalParameters.names() = startingValues_.names();
    
    if(!fitResults_.convergence) Rcpp::warning("Optimizer did not converge");
    Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("fit") = fitResults_.fit,
      Rcpp::Named("convergence") = fitResults_.convergence,
      Rcpp::Named("rawParameters") = finalParameters,
      Rcpp::Named("fits") = fitResults_.fits
    );
    return(result);
  }
};

RCPP_EXPOSED_CLASS(istaLspGeneralPurpose)
RCPP_MODULE(istaLspGeneralPurpose_cpp){
  using namespace Rcpp;
  Rcpp::class_<istaLspGeneralPurpose>( "istaLspGeneralPurpose" )
  .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new istaLspGeneralPurpose.")
  // methods
  .method( "optimize", &istaLspGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values, theta, and lambda")
  ;
}

RCPP_EXPOSED_CLASS(istaLspGeneralPurposeCpp)
  RCPP_MODULE(istaLspGeneralPurposeCpp_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaLspGeneralPurposeCpp>( "istaLspGeneralPurposeCpp" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new istaLspGeneralPurposeCpp.")
    // methods
    .method( "optimize", &istaLspGeneralPurposeCpp::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values, theta, and lambda")
    ;
  }