#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "gpFitFramework.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class istaEnetGeneralPurpose{
public:
  
  
  // settings
  Rcpp::NumericVector startingValues;
  double lambda, alpha;
  const Rcpp::NumericVector weights;
  
  // control
  const int i_k;
  const double L0;
  const double eta;
  const int maxIterOut;
  const int maxIterIn;
  const double breakOuter;
  
  // constructor
  istaEnetGeneralPurpose(
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
      double lambda_,
      double alpha_){
    
    generalPurposeFitFramework gpFF(fitFunction, gradientFunction, userSuppliedElements);
    
    startingValues = startingValues_;
    lambda = lambda_;
    alpha = alpha_;
    
    linr::tuningParametersEnet tp;
    tp.lambda = lambda;
    tp.weights = weights;
    
    linr::proximalOperatorLasso proximalOperatorLasso_;
    linr::penaltyLASSO penalty_;
    linr::penaltyRidge smoothPenalty_;
    
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
      smoothPenalty_,
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

RCPP_EXPOSED_CLASS(istaEnetGeneralPurpose)
  RCPP_MODULE(istaEnetGeneralPurpose_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaEnetGeneralPurpose>( "istaEnetGeneralPurpose" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new istaEnetGeneralPurpose.")
      .field_readonly( "lambda", &istaEnetGeneralPurpose::lambda, "tuning parameter lambda")
    // methods
    .method( "optimize", &istaEnetGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values and lambda")
    ;
  }

