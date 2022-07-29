#include <RcppArmadillo.h>
#include "SEM.hpp"
#include "SEMFitFramework.hpp"

// [[Rcpp :: depends ( RcppArmadillo )]]

//'@name istaEnet
//'@title elastic net optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a lambda and an alpha value.
//'@returns a list with fit results
class istaEnet{
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
  istaEnet(
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
      double lambda_, 
      double alpha_){
    
    SEMFitFramework SEMFF(SEM_);
    
    int sampleSize = SEMFF.SEM.rawData.n_rows;
    
    lessSEM::tuningParametersEnet tp;
    tp.alpha = alpha_;
    tp.lambda = lambda_;
    tp.weights = weights;
    lessSEM::tuningParametersEnet smoothTp = tp;
    
    lessSEM::proximalOperatorLasso proximalOperatorLasso_;
    lessSEM::penaltyLASSO penalty_;
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
      proximalOperatorLasso_,
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

//'@name glmnetEnet
//'@title elastic net optimization with glmnet optimizer
//'@description Object for elastic net optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a lambda and an alpha value.
//'@returns a list with fit results
//'
class glmnetEnet{
public:
  
  Rcpp::NumericVector startingValues;
  double alpha;
  double lambda;
  const arma::rowvec weights;
  // control optimizer
  arma::mat initialHessian;
  const double stepSize;
  const double sigma;
  const double gamma;
  const int maxIterOut; // maximal number of outer iterations
  const int maxIterIn; // maximal number of inner iterations
  const int maxIterLine;
  const double breakOuter; // change in fit required to break the outer iteration
  const double breakInner;
  const lessSEM::convergenceCriteriaGlmnet convergenceCriterion; // this is related to the inner
  // breaking condition. 
  const int verbose;
  
  
  // constructor
  glmnetEnet(
    const arma::rowvec weights_,
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
    
    convergenceCriterion(static_cast<lessSEM::convergenceCriteriaGlmnet>(Rcpp::as<int> (control["convergenceCriterion"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  void setHessian(arma::mat newHessian){
    initialHessian = newHessian;
  }
  
  Rcpp::List optimize(Rcpp::NumericVector startingValues_, 
                      SEMCpp& SEM_, 
                      double lambda_, 
                      double alpha_){
    
    SEMFitFramework SEMFF(SEM_);
    
    int N = SEMFF.SEM.rawData.n_rows;
    
    lessSEM::tuningParametersEnet tp;
    tp.alpha = alpha_;
    tp.lambda = lambda_*N;
    tp.weights = weights;
    
    lessSEM::proximalOperatorLasso proximalOperatorLasso_;
    lessSEM::penaltyLASSO penalty_;
    lessSEM::penaltyRidge smoothPenalty_;
    
    lessSEM::controlGLMNET control_ = {
      initialHessian,
      stepSize,
      sigma,
      gamma,
      maxIterOut, // maximal number of outer iterations
      maxIterIn, // maximal number of inner iterations
      maxIterLine,
      breakOuter*N, // change in fit required to break the outer iteration
      breakInner*N,
      convergenceCriterion, // this is related to the inner
      // breaking condition. 
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::glmnet(
      SEMFF,
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

//'@name bfgsEnet
//'@title smoothly approximated elastic net
//'@description Object for smoothly approximated elastic net optimization with
//'bfgs optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a lambda and an alpha value.
//'@returns a list with fit results
//'
class bfgsEnet{
public:
  
  Rcpp::NumericVector startingValues;
  double alpha;
  double lambda;
  const arma::rowvec weights;
  // control optimizer
  double epsilon; // controls smoothness of lasso penalty function
  arma::mat initialHessian;
  const double stepSize;
  const double sigma;
  const double gamma;
  const int maxIterOut; // maximal number of outer iterations
  const int maxIterIn; // maximal number of inner iterations
  const int maxIterLine;
  const double breakOuter; // change in fit required to break the outer iteration
  const double breakInner;
  const lessSEM::convergenceCriteriaBFGS convergenceCriterion; // this is related to the inner
  // breaking condition. 
  const int verbose;
  
  
  // constructor
  bfgsEnet(
    const arma::rowvec weights_,
    Rcpp::List control
  ): 
    weights(weights_),
    epsilon(Rcpp::as<double> (control["epsilon"])),
    initialHessian(Rcpp::as<arma::mat> (control["initialHessian"])),
    stepSize(Rcpp::as<double> (control["stepSize"])),
    sigma(Rcpp::as<double> (control["sigma"])),
    gamma(Rcpp::as<double> (control["gamma"])),
    
    maxIterOut(Rcpp::as<int> (control["maxIterOut"])),
    maxIterIn(Rcpp::as<int> (control["maxIterIn"])),
    maxIterLine(Rcpp::as<int> (control["maxIterLine"])),
    
    breakOuter(Rcpp::as<double> (control["breakOuter"])),
    breakInner(Rcpp::as<double> (control["breakInner"])),
    
    convergenceCriterion(static_cast<lessSEM::convergenceCriteriaBFGS>(Rcpp::as<int> (control["convergenceCriterion"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  void setHessian(arma::mat newHessian){
    initialHessian = newHessian;
  }
  
  Rcpp::List optimize(Rcpp::NumericVector startingValues_, 
                      SEMCpp& SEM_, 
                      double lambda_, 
                      double alpha_){
    
    SEMFitFramework SEMFF(SEM_);
    
    int N = SEMFF.SEM.rawData.n_rows;
    
    lessSEM::tuningParametersSmoothElasticNet tp;
    tp.epsilon = epsilon;
    tp.alpha = alpha_;
    tp.lambda = lambda_*N;
    tp.weights = weights;
    
    lessSEM::smoothElasticNet smoothPenalty_;
    
    lessSEM::controlBFGS control_ = {
      initialHessian,
      stepSize,
      sigma,
      gamma,
      maxIterOut, // maximal number of outer iterations
      maxIterIn, // maximal number of inner iterations
      maxIterLine,
      breakOuter*N, // change in fit required to break the outer iteration
      breakInner*N,
      convergenceCriterion, // this is related to the inner
      // breaking condition. 
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::bfgsOptim(
      SEMFF,
      startingValues_,
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

RCPP_EXPOSED_CLASS(bfgsEnet)
  RCPP_MODULE(bfgsEnet_cpp){
    using namespace Rcpp;
    Rcpp::class_<bfgsEnet>( "bfgsEnet" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaEnet.")
    // methods
    .method( "setHessian", &bfgsEnet::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &bfgsEnet::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

RCPP_EXPOSED_CLASS(glmnetEnet)
  RCPP_MODULE(glmnetEnet_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetEnet>( "glmnetEnet" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaEnet.")
    // methods
    .method( "setHessian", &glmnetEnet::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &glmnetEnet::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

RCPP_EXPOSED_CLASS(istaEnet)
  RCPP_MODULE(istaEnet_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaEnet>( "istaEnet" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaEnet.")
    // methods
    .method( "optimize", &istaEnet::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

