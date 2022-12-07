#include <RcppArmadillo.h>
#include "SEM.h"
#include "mgSEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

//' mcpPenalty_C
//' 
//' @param par single parameter value
//' @param lambda_p lambda value for this parameter
//' @param theta theta value for this parameter
//' @returns penalty value
// [[Rcpp::export]]
double mcpPenalty_C(const double par,
                     const double lambda_p,
                     const double theta){
  return(lessSEM::mcpPenalty(par, lambda_p, theta));
}

template<typename sem>
class istaMcp{
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
  istaMcp(
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
    
    lessSEM::tuningParametersMcp tp;
    tp.theta = theta_;
    tp.lambda = lambda_;
    tp.weights = weights;
    
    // we won't need the smooth penalty; but we need to specify some tuning 
    // parameters for the function call
    lessSEM::tuningParametersMcp smoothTp;
    smoothTp.lambda = 0.0;
    
    lessSEM::proximalOperatorMcp proximalOperatorMcp_;
    lessSEM::penaltyMcp penalty_;
    lessSEM::noSmoothPenalty<lessSEM::tuningParametersMcp> smoothPenalty_;
    
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
      proximalOperatorMcp_,
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

//'@name istaMcpSEM
//'@title mcp optimization with ista
//'@description Object for mcp optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
typedef istaMcp<SEMCpp> istaMcpSEM;
RCPP_EXPOSED_CLASS_NODECL(istaMcpSEM)
RCPP_MODULE(istaMcpSEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<istaMcpSEM>( "istaMcpSEM" )
  .constructor<arma::rowvec,Rcpp::List>("Creates a new istaMcpSEM.")
  // methods
  .method( "optimize", &istaMcpSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
  ;
}

//'@name istaMcpMgSEM
//'@title mcp optimization with ista
//'@description Object for mcp optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
typedef istaMcp<mgSEM> istaMcpMgSEM;
RCPP_EXPOSED_CLASS_NODECL(istaMcpMgSEM)
  RCPP_MODULE(istaMcpMgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaMcpMgSEM>( "istaMcpMgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaMcpMgSEM.")
    // methods
    .method( "optimize", &istaMcpMgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }