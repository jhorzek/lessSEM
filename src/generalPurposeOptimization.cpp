#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "gpFitFramework.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

class istaEnetGeneralPurpose{
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
  const linr::convCritInnerIsta convCritInner;
  const double sigma;
  const linr::stepSizeInheritance stepSizeInh;
  const int verbose; 
  
  // constructor
  istaEnetGeneralPurpose(
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
    convCritInner(static_cast<linr::convCritInnerIsta>(Rcpp::as<int> (control["convCritInner"]))),
    sigma(Rcpp::as<double> (control["sigma"])),
    stepSizeInh(static_cast<linr::stepSizeInheritance>(Rcpp::as<int> (control["stepSizeInheritance"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      Rcpp::Function fitFunction,
      Rcpp::Function gradientFunction,
      Rcpp::List userSuppliedElements,
      double lambda_,
      double alpha_){
    
    generalPurposeFitFramework gpFF(fitFunction, gradientFunction, userSuppliedElements);
    
    linr::tuningParametersEnet tp;
    tp.lambda = lambda_;
    tp.weights = weights;
    tp.alpha = alpha_;
    
    linr::proximalOperatorLasso proximalOperatorLasso_;
    linr::penaltyLASSO penalty_;
    linr::penaltyRidge smoothPenalty_;
    
    linr::control controlIsta = {
      L0,
      eta,
      accelerate,
      maxIterOut,
      maxIterIn,
      breakOuter,
      convCritInner,
      sigma,
      stepSizeInh,
      verbose
    };
    
    linr::fitResults fitResults_ = linr::ista(
      gpFF,
      startingValues_,
      proximalOperatorLasso_,
      penalty_,
      smoothPenalty_,
      tp,
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

class glmnetEnetGeneralPurpose{
public:
  
  
  Rcpp::NumericVector startingValues;
  double alpha;
  double lambda;
  const arma::rowvec weights;
  // control optimizer
  const arma::mat initialHessian;
  const double stepSize;
  const double sigma;
  const double gamma;
  const int maxIterOut; // maximal number of outer iterations
  const int maxIterIn; // maximal number of inner iterations
  const int maxIterLine;
  const double breakOuter; // change in fit required to break the outer iteration
  const double breakInner;
  const linr::convergenceCriteriaGlmnet convergenceCriterion; // this is related to the inner
  // breaking condition. 
  const int verbose;
  
  // constructor
  glmnetEnetGeneralPurpose(
    const Rcpp::NumericVector weights_,
    Rcpp::List control
  ):
    weights(weights_),
    initialHessian(Rcpp::as<arma::mat> (control["initialHessian"])),
    stepSize(Rcpp::as<double> (control["stepSize"])),
    sigma(Rcpp::as<double> (control["sigma"])),
    gamma(Rcpp::as<double> (control["gamma"])),
    
    maxIterOut(Rcpp::as<int> (control["maxIterOut"])),
    maxIterIn(Rcpp::as<int> (control["maxIterIn"])),
    maxIterLine(Rcpp::as<int> (control["maxIterLine"])),
    
    breakOuter(Rcpp::as<double> (control["breakOuter"])),
    breakInner(Rcpp::as<double> (control["breakInner"])),
    
    convergenceCriterion(static_cast<linr::convergenceCriteriaGlmnet>(Rcpp::as<int> (control["convergenceCriterion"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      Rcpp::Function fitFunction,
      Rcpp::Function gradientFunction,
      Rcpp::List userSuppliedElements,
      double lambda_,
      double alpha_){
    
    generalPurposeFitFramework gpFF(fitFunction, 
                                    gradientFunction,
                                    userSuppliedElements);
    
    linr::tuningParametersEnet tp;
    tp.lambda = lambda_;
    tp.weights = weights;
    tp.alpha = alpha_;
    
    linr::proximalOperatorLasso proximalOperatorLasso_;
    linr::penaltyLASSO penalty_;
    linr::penaltyRidge smoothPenalty_;
    
    linr::controlGLMNET control_ = {
      initialHessian,
      stepSize,
      sigma,
      gamma,
      maxIterOut, // maximal number of outer iterations
      maxIterIn, // maximal number of inner iterations
      maxIterLine,
      breakOuter, // change in fit required to break the outer iteration
      breakInner,
      convergenceCriterion, // this is related to the inner
      // breaking condition. 
      verbose
    };
    
    linr::fitResults fitResults_ = linr::glmnet(
      gpFF,
      startingValues_,
      penalty_,
      smoothPenalty_,
      tp,
      control_
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
      Rcpp::Named("fits") = fitResults_.fits,
      Rcpp::Named("Hessian") = fitResults_.Hessian
    );
    return(result);
  }
};

RCPP_EXPOSED_CLASS(glmnetEnetGeneralPurpose)
  RCPP_MODULE(glmnetEnetGeneralPurpose_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetEnetGeneralPurpose>( "glmnetEnetGeneralPurpose" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new glmnetEnetGeneralPurpose.")
    // methods
    .method( "optimize", &glmnetEnetGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values and lambda")
    ;
  }

RCPP_EXPOSED_CLASS(istaEnetGeneralPurpose)
  RCPP_MODULE(istaEnetGeneralPurpose_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaEnetGeneralPurpose>( "istaEnetGeneralPurpose" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new istaEnetGeneralPurpose.")
    // methods
    .method( "optimize", &istaEnetGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values and lambda")
    ;
  }

