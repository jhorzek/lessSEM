#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "SEMFitFramework.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class istaCappedL1{
  public:
    
    Rcpp::NumericVector startingValues;
  double alpha;
  double lambda;
  const arma::rowvec weights;
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
  istaCappedL1(
    const arma::rowvec weights_,
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
    SEMCpp& SEM_,
    double theta_,
    double lambda_, 
    double alpha_){
    
    SEMFitFramework SEMFF(SEM_);
    
    int N = SEMFF.SEM.rawData.n_rows;
    
    lessSEM::tuningParametersCappedL1 tp;
    tp.theta = theta_;
    tp.alpha = alpha_;
    tp.lambda = lambda_*N;
    tp.weights = weights;
    
    lessSEM::tuningParametersEnet smoothTp;
    smoothTp.alpha = alpha_;
    smoothTp.lambda = lambda_*N;
    smoothTp.weights = weights;
    
    lessSEM::proximalOperatorCappedL1 proximalOperatorCappedL1_;
    lessSEM::penaltyCappedL1 penalty_;
    lessSEM::penaltyRidge smoothPenalty_;
    
    lessSEM::control controlIsta = {
      L0,
      eta,
      accelerate,
      maxIterOut,
      maxIterIn,
      breakOuter*N,
      convCritInner,
      sigma,
      stepSizeInh,
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::ista(
      SEMFF,
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

RCPP_EXPOSED_CLASS(istaCappedL1)
  RCPP_MODULE(istaCappedL1_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaCappedL1>( "istaCappedL1" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaCappedL1.")
    // methods
    .method( "optimize", &istaCappedL1::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }