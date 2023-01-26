#include <RcppArmadillo.h>
#include "mgSEM.h"
#include "SEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]


template<typename sem>
class istaMixedPenalty{
public:
  
  Rcpp::NumericVector startingValues;
  
  std::vector<lessSEM::penaltyType> pType;
  arma::rowvec lambda;
  arma::rowvec theta;
  arma::rowvec alpha;
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
  istaMixedPenalty(
    const arma::rowvec weights_,
    const std::vector<int> pType_int,
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
    verbose(Rcpp::as<int> (control["verbose"])){
    
    pType.resize(pType_int.size());
    
    for(int i = 0; i < pType_int.size(); i++){
      pType.at(i) = static_cast<lessSEM::penaltyType>(pType_int.at(i));
    }
    
  }
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      sem& SEM_,
      arma::rowvec theta_,
      arma::rowvec lambda_, 
      arma::rowvec alpha_){
    
    SEMFitFramework<sem> SEMFF(SEM_);
    
    int sampleSize = SEMFF.SEM.sampleSize;
    
    lessSEM::tuningParametersMixedPenalty tp;
    tp.pt = pType;
    tp.lambda = lambda_;
    tp.theta = theta_;
    tp.alpha = alpha_;
    tp.weights = weights;
    
    lessSEM::tuningParametersEnet smoothTp;
    smoothTp.alpha = 0.0;
    smoothTp.lambda = 0.0;
    smoothTp.weights = weights;
    
    lessSEM::proximalOperatorMixedPenalty proximalOperatorMixedPenalty_;
    lessSEM::penaltyMixedPenalty penalty_;
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
      proximalOperatorMixedPenalty_,
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

//'@name istaMixedPenaltySEM
//'@title mixed penalty optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta value, a lambda and an alpha value (alpha must be 1).
//'@returns a list with fit results
typedef istaMixedPenalty<SEMCpp> istaMixedPenaltySEM;
RCPP_EXPOSED_CLASS_NODECL(istaMixedPenaltySEM)
  RCPP_MODULE(istaMixedPenaltySEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaMixedPenaltySEM>( "istaMixedPenaltySEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaMixedPenaltySEM.")
    // methods
    .method( "optimize", &istaMixedPenaltySEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }

//'@name istaMixedPenaltymgSEM
//'@title mixed penalty optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta value, a lambda and an alpha value (alpha must be 1).
//'@returns a list with fit results
typedef istaMixedPenalty<mgSEM> istaMixedPenaltymgSEM;
RCPP_EXPOSED_CLASS_NODECL(istaMixedPenaltymgSEM)
  RCPP_MODULE(istaMixedPenaltymgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaMixedPenaltymgSEM>( "istaMixedPenaltymgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaMixedPenaltymgSEM.")
    // methods
    .method( "optimize", &istaMixedPenaltymgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }