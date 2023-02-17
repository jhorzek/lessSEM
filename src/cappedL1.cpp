#include <RcppArmadillo.h>
#include "mgSEM.h"
#include "SEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]


template<typename sem>
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
      sem& SEM_,
      double theta_,
      double lambda_, 
      double alpha_){
    
    if(alpha_ != 1.0) Rcpp::stop("alpha must be set to 1.");
    
    SEMFitFramework<sem> SEMFF(SEM_);
    
    int sampleSize = SEMFF.SEM.sampleSize;
    
    lessSEM::tuningParametersCappedL1 tp;
    tp.theta = theta_;
    tp.alpha = alpha_;
    tp.lambda = lambda_;
    tp.weights = weights;
    
    lessSEM::tuningParametersEnet smoothTp;
    smoothTp.alpha = alpha_;
    smoothTp.lambda = lambda_;
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
      proximalOperatorCappedL1_,
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

//'@name istaCappedL1SEM
//'@title cappedL1 optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta value, a lambda and an alpha value (alpha must be 1).
//'@returns a list with fit results
typedef istaCappedL1<SEMCpp> istaCappedL1SEM;
RCPP_EXPOSED_CLASS_NODECL(istaCappedL1SEM)
  RCPP_MODULE(istaCappedL1SEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaCappedL1SEM>( "istaCappedL1SEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaCappedL1SEM.")
    // methods
    .method( "optimize", &istaCappedL1SEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }

//'@name istaCappedL1mgSEM
//'@title cappedL1 optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta value, a lambda and an alpha value (alpha must be 1).
//'@returns a list with fit results
typedef istaCappedL1<mgSEM> istaCappedL1mgSEM;
RCPP_EXPOSED_CLASS_NODECL(istaCappedL1mgSEM)
  RCPP_MODULE(istaCappedL1mgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaCappedL1mgSEM>( "istaCappedL1mgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaCappedL1mgSEM.")
    // methods
    .method( "optimize", &istaCappedL1mgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }