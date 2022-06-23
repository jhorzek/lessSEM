#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "gpFitFramework.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class istaLASSOGeneralPurpose{
public:
  
  
  // settings
  Rcpp::NumericVector startingValues;
  double lambda;
  const Rcpp::NumericVector weights;
  
  // control
  const int i_k;
  const double L0;
  const double eta;
  const int maxIterOut;
  const int maxIterIn;
  const double breakOuter;
  
  // constructor
  istaLASSOGeneralPurpose(
    const Rcpp::NumericVector weights_,
    Rcpp::List control
  ):
    weights(weights_),
    i_k(Rcpp::as<int>(control["i_k"])),
    L0(Rcpp::as<double> (control["L0"])),
    eta(Rcpp::as<double> (control["eta"])),
    maxIterOut(Rcpp::as<int> (control["maxIterOut"])),
    maxIterIn(Rcpp::as<int> (control["maxIterIn"])),
    breakOuter(Rcpp::as<double> (control["breakOuter"])){}
  
  Rcpp::List optimize(
      Rcpp::Function fitFunction,
      Rcpp::Function gradientFunction,
      Rcpp::List userSuppliedElements,
      Rcpp::NumericVector startingValues_, 
      double lambda_){
    
    generalPurposeFitFramework gpFF(fitFunction, gradientFunction, userSuppliedElements);
    
    startingValues = startingValues_;
    lambda = lambda_;
    
    linr::tuningParametersLasso tp;
    tp.lambda = lambda;
    tp.weights = weights;
    
    linr::proximalOperatorLasso proximalOperatorLasso_;
    linr::penaltyLASSO penalty_;
    
    linr::control controlIsta = {
      i_k,
      L0,
      eta,
      maxIterOut,
      maxIterIn,
      breakOuter
    };
    
    linr::fitResults fitResults_ = linr::ista(
      gpFF,
      startingValues,
      proximalOperatorLasso_,
      penalty_,
      tp,
      controlIsta
    );
    
    if(!fitResults_.convergence) Rcpp::warning("Optimizer did not converge");
    Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("fit") = fitResults_.fit,
      Rcpp::Named("convergence") = fitResults_.convergence,
      Rcpp::Named("rawParameters") = fitResults_.parameterValues,
      Rcpp::Named("fits") = fitResults_.fits
    );
    return(result);
  }
};

RCPP_EXPOSED_CLASS(istaLASSOGeneralPurpose)
  RCPP_MODULE(istaLASSOGeneralPurpose_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaLASSOGeneralPurpose>( "istaLASSOGeneralPurpose" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new istaLASSOGeneralPurpose.")
      .field_readonly( "lambda", &istaLASSOGeneralPurpose::lambda, "tuning parameter lambda")
    // methods
    .method( "optimize", &istaLASSOGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values and lambda")
    ;
  }

