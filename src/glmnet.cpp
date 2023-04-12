#include <RcppArmadillo.h>
#include "SEM.h"
#include "mgSEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

template<typename sem>
class glmnetScad{
public:
  
  Rcpp::NumericVector startingValues;
  arma::rowvec weights;
  double theta;
  double lambda;
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
  glmnetScad(
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
    verbose(Rcpp::as<int> (control["verbose"])){
    
    for(auto w: weights){
      if((w != 0.0) && (w != 1.0))
        Rcpp::stop("All weights must be either 0 or 1");
    }
  }
  
  void setHessian(arma::mat newHessian){
    initialHessian = newHessian;
  }
  
  Rcpp::List optimize(Rcpp::NumericVector startingValues_, 
                      sem& SEM_, 
                      double theta_,
                      double lambda_){
    
    const double N = SEM_.sampleSize;
    
    SEMFitFramework<sem> SEMFF(SEM_, 1.0/N);
    
    lessSEM::tuningParametersScadGlmnet tp;
    
    tp.lambda = lambda_;
    tp.theta = theta_;
    tp.weights = weights;
    
    lessSEM::penaltySCADGlmnet penalty_;
    lessSEM::noSmoothPenalty<lessSEM::tuningParametersScadGlmnet> smoothPenalty_;
    
    lessSEM::controlGLMNET control_ = {
      initialHessian/N,
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
    
    lessSEM::fitResults fitResults_ = lessSEM::glmnet(
      SEMFF,
      startingValues_,
      penalty_,
      smoothPenalty_,
      tp,
      control_
    );
    
    Rcpp::NumericVector finalParameters(fitResults_.parameterValues.n_elem);
    for(unsigned int i = 0; i < fitResults_.parameterValues.n_elem; i++){
      finalParameters.at(i) = fitResults_.parameterValues.at(i);
    }
    finalParameters.names() = startingValues_.names();
    
    if(!fitResults_.convergence) Rcpp::warning("Optimizer did not converge");
    Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("fit") = N*fitResults_.fit,
      Rcpp::Named("convergence") = fitResults_.convergence,
      Rcpp::Named("rawParameters") = finalParameters,
      Rcpp::Named("fits") = N*fitResults_.fits,
      Rcpp::Named("Hessian") = N*fitResults_.Hessian
    );
    return(result);
  }
};

//'@name glmnetScadSEM
//'@title scad optimization with glmnet optimizer
//'@description Object for scad optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef glmnetScad<SEMCpp> glmnetScadSEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetScadSEM)
  RCPP_MODULE(glmnetScadSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetScadSEM>( "glmnetScadSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new glmnetScadSEM.")
    // methods
    .method( "setHessian", &glmnetScadSEM::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &glmnetScadSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

//'@name glmnetScadMgSEM
//'@title scad optimization with glmnet optimizer
//'@description Object for scad optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef glmnetScad<mgSEM> glmnetScadMgSEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetScadMgSEM)
  RCPP_MODULE(glmnetScadMgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetScadMgSEM>( "glmnetScadMgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new glmnetScadMgSEM.")
    // methods
    .method( "setHessian", &glmnetScadMgSEM::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &glmnetScadMgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }
