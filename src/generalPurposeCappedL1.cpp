#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "gpFitFramework.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class cappedL1EnetGeneralPurpose{
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
  cappedL1EnetGeneralPurpose(
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
    
    lessSEM::tuningParametersCappedL1 tp;
    tp.theta = theta_;
    tp.lambda = lambda_;
    tp.weights = weights;
    
    // we won't need the smooth penalty; but we need to specify some tuning 
    // parameters for the function call
    lessSEM::tuningParametersCappedL1 smoothTp;
    smoothTp.lambda = 0.0;
    
    lessSEM::proximalOperatorCappedL1 proximalOperatorCappedL1_;
    lessSEM::penaltyLSP penalty_;
    lessSEM::noSmoothPenalty<lessSEM::tuningParametersCappedL1> smoothPenalty_;
    
    int sampleSize = 1;
    
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
      proximalOperatorCappedL1_,
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

RCPP_EXPOSED_CLASS(cappedL1EnetGeneralPurpose)
  RCPP_MODULE(cappedL1EnetGeneralPurpose_cpp){
    using namespace Rcpp;
    Rcpp::class_<cappedL1EnetGeneralPurpose>( "cappedL1EnetGeneralPurpose" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new cappedL1EnetGeneralPurpose.")
    // methods
    .method( "optimize", &cappedL1EnetGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values, theta, and lambda")
    ;
  }