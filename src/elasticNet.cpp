#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "SEMFitFramework.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class istaEnet{
public:
  
  Rcpp::NumericVector startingValues;
  double alpha;
  double lambda;
  const Rcpp::NumericVector weights;
  // control optimizer
  const double L0;
  const double eta;
  const int maxIterOut;
  const int maxIterIn;
  const double breakOuter;
  const linr::convCritInnerIsta convCritInner;
  const double sigma;
  const linr::stepSizeInheritance stepSizeInh;
  const int verbose; 
  
  
  // constructor
  istaEnet(
    const Rcpp::NumericVector weights_,
    Rcpp::List control
  ): 
    weights(weights_),
    L0(Rcpp::as<double> (control["L0"])),
    eta(Rcpp::as<double> (control["eta"])),
    maxIterOut(Rcpp::as<int> (control["maxIterOut"])),
    maxIterIn(Rcpp::as<int> (control["maxIterIn"])),
    breakOuter(Rcpp::as<double> (control["breakOuter"])),
    convCritInner(static_cast<linr::convCritInnerIsta>(Rcpp::as<int> (control["convCritInner"]))),
    sigma(Rcpp::as<double> (control["sigma"])),
    stepSizeInh(static_cast<linr::stepSizeInheritance>(Rcpp::as<int> (control["stepSizeInheritance"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  Rcpp::List optimize(SEMCpp& SEM_, 
                      Rcpp::NumericVector startingValues_, 
                      double lambda_, 
                      double alpha_){
    
    SEMFitFramework SEMFF(SEM_);
    
    int N = SEMFF.SEM.rawData.n_rows;
    
    linr::tuningParametersEnet tp;
    tp.alpha = alpha_;
    tp.lambda = lambda_*N;
    tp.weights = weights;
    
    linr::proximalOperatorLasso proximalOperatorLasso_;
    linr::penaltyLASSO penalty_;
    linr::penaltyRidge smoothPenalty_;
    
    linr::control controlIsta = {
      L0,
      eta,
      maxIterOut,
      maxIterIn,
      breakOuter*N,
      convCritInner,
      sigma,
      stepSizeInh,
      verbose
    };
    
    linr::fitResults fitResults_ = linr::ista(
      SEMFF,
      startingValues_,
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

RCPP_EXPOSED_CLASS(istaEnet)
  RCPP_MODULE(istaEnet_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaEnet>( "istaEnet" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new istaEnet.")
      .field_readonly( "lambda", &istaEnet::lambda, "tuning parameter lambda")
    // methods
    .method( "optimize", &istaEnet::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values and lambda")
    ;
  }

