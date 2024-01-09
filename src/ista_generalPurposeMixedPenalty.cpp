#include <RcppArmadillo.h>
#include "SEM.h"
#include "gpFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

class istaMixedPenaltyGeneralPurpose{
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
  istaMixedPenaltyGeneralPurpose(
    const arma::rowvec weights_,
    const std::vector<std::string> pType_str,
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
    
    pType = lessSEM::stringPenaltyToPenaltyType(pType_str);
    
  }
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      Rcpp::Function fitFunction,
      Rcpp::Function gradientFunction,
      Rcpp::List userSuppliedElements,
    arma::rowvec lambda_, 
    arma::rowvec theta_,
    arma::rowvec alpha_){
    
    generalPurposeFitFramework gpFF(fitFunction, 
                                    gradientFunction, 
                                    userSuppliedElements);
    
    int sampleSize = 1;
    
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
    
    initializeMixedProximalOperators(proximalOperatorMixedPenalty_,
                                     pType);
    initializeMixedPenalties(penalty_, 
                             pType);
    
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
      proximalOperatorMixedPenalty_,
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

class istaMixedPenaltyGeneralPurposeCpp{
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
  istaMixedPenaltyGeneralPurposeCpp(
    const arma::rowvec weights_,
    const std::vector<std::string> pType_str,
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
    
    pType = lessSEM::stringPenaltyToPenaltyType(pType_str);
    
  }
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      SEXP fitFunction,
      SEXP gradientFunction,
      Rcpp::List userSuppliedElements,
      arma::rowvec lambda_, 
      arma::rowvec theta_,
      arma::rowvec alpha_){
    
    generalPurposeFitFrameworkCpp gpFF(startingValues_, 
                                       fitFunction, 
                                       gradientFunction, 
                                       userSuppliedElements);
    
    int sampleSize = 1;
    
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
    
    initializeMixedProximalOperators(proximalOperatorMixedPenalty_,
                                     pType);
    initializeMixedPenalties(penalty_, 
                             pType);
    
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
      proximalOperatorMixedPenalty_,
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

//'@name istaMixedPenaltyGeneralPurpose
//'@title mixed penalty optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. 
//'@field optimize optimize the model. 
//'@returns a list with fit results
RCPP_EXPOSED_CLASS_NODECL(istaMixedPenaltyGeneralPurpose)
  RCPP_MODULE(istaMixedPenaltyGeneralPurpose_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaMixedPenaltyGeneralPurpose>( "istaMixedPenaltyGeneralPurpose" )
      .constructor<arma::rowvec,std::vector<std::string>,Rcpp::List>("Creates a new istaMixedPenaltyGeneralPurpose.")
    // methods
    .method( "optimize", &istaMixedPenaltyGeneralPurpose::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
    ;
  }

//'@name istaMixedPenaltyGeneralPurposeCpp
//'@title mixed penalty optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter, (2) a vector indicating which penalty is used, and (3) a list with control elements
//'@field optimize optimize the model. 
//'@returns a list with fit results
RCPP_EXPOSED_CLASS_NODECL(istaMixedPenaltyGeneralPurposeCpp)
RCPP_MODULE(istaMixedPenaltyGeneralPurposeCpp_cpp){
  using namespace Rcpp;
  Rcpp::class_<istaMixedPenaltyGeneralPurposeCpp>( "istaMixedPenaltyGeneralPurposeCpp" )
  .constructor<arma::rowvec,std::vector<std::string>,Rcpp::List>("Creates a new istaMixedPenaltyGeneralPurposeCpp.")
  // methods
  .method( "optimize", &istaMixedPenaltyGeneralPurposeCpp::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, theta, lambda, and alpha")
  ;
}