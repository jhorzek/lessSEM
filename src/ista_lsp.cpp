#include <RcppArmadillo.h>
#include "SEM.h"
#include "mgSEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

template<typename sem>
class istaLSP{
public:
  
  Rcpp::NumericVector startingValues;
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
  istaLSP(
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
      sem& SEM_,
      double theta_,
      double lambda_){
    
    SEMFitFramework<sem> SEMFF(SEM_);
    
    int sampleSize = SEMFF.SEM.sampleSize;
    
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
      SEMFF,
      startingValues_,
      proximalOperatorLSP_,
      penalty_,
      smoothPenalty_,
      tp,
      smoothTp,
      controlIsta
    );
    
    Rcpp::NumericVector finalParameters(fitResults_.parameterValues.n_elem);
    for(unsigned int i = 0; i < fitResults_.parameterValues.n_elem; i++){
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

//'@name istaLSPSEM
//'@title lsp optimization with ista
//'@description Object for lsp optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef istaLSP<SEMCpp> istaLSPSEM;
RCPP_EXPOSED_CLASS_NODECL(istaLSPSEM)
  RCPP_MODULE(istaLSPSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaLSPSEM>( "istaLSPSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaLSPSEM.")
    // methods
    .method( "optimize", &istaLSPSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }

//'@name istaLSPMgSEM
//'@title lsp optimization with ista
//'@description Object for lsp optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef istaLSP<mgSEM> istaLSPMgSEM;
RCPP_EXPOSED_CLASS_NODECL(istaLSPMgSEM)
  RCPP_MODULE(istaLSPMgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaLSPMgSEM>( "istaLSPMgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaLSPMgSEM.")
    // methods
    .method( "optimize", &istaLSPMgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }