#include <RcppArmadillo.h>
#include "SEM.h"
#include "mgSEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

template<typename sem>
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
      sem& SEM_,
      double lambda_, 
      double alpha_){
    
    SEMFitFramework<sem> SEMFF(SEM_);
    
    int sampleSize = SEMFF.SEM.sampleSize;
    
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

template<typename sem>
class glmnetEnet{
public:
  
  Rcpp::NumericVector startingValues;
  arma::rowvec alpha;
  arma::rowvec lambda;
  arma::rowvec weights;
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
                      sem& SEM_, 
                      arma::rowvec lambda_, 
                      arma::rowvec alpha_){
    
    SEMFitFramework<sem> SEMFF(SEM_);
    
    int N = SEMFF.SEM.sampleSize;
    
    lessSEM::tuningParametersEnetGlmnet tp;
    
    tp.weights = weights;
    
    if((alpha_.n_elem != tp.weights.n_elem) && (alpha_.n_elem == 1)){
      // glmnet wants one alpha for each parameter; assume that they are all the
      // same value if the user just supplied one.
      tp.alpha = arma::rowvec(tp.weights.n_elem);
      tp.alpha.fill(alpha_.at(0));
    }else if(alpha_.n_elem == tp.weights.n_elem){
      tp.alpha = alpha_;
    }else{
      Rcpp::stop("alpha must be either of size 1 or of the same length as the weights.");
    }
    
    if((lambda_.n_elem != tp.weights.n_elem) && (lambda_.n_elem == 1)){
      // glmnet wants one alpha for each parameter; assume that they are all the
      // same value if the user just supplied one.
      tp.lambda = arma::rowvec(tp.weights.n_elem);
      tp.lambda.fill(lambda_.at(0)*N);
    }else if(lambda_.n_elem == tp.weights.n_elem){
      tp.lambda = lambda_*N;
    }else{
      Rcpp::stop("lambda must be either of size 1 or of the same length as the weights.");
    }
    
    lessSEM::penaltyLASSOGlmnet penalty_;
    lessSEM::penaltyRidgeGlmnet smoothPenalty_;
    
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
    for(unsigned int i = 0; i < fitResults_.parameterValues.n_elem; i++){
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
template<typename sem>
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
                      sem& SEM_, 
                      double lambda_, 
                      double alpha_){
    
    SEMFitFramework<sem> SEMFF(SEM_);
    
    int N = SEMFF.SEM.sampleSize;
    
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
    for(unsigned int i = 0; i < fitResults_.parameterValues.n_elem; i++){
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

//'@name bfgsEnetSEM
//'@title smoothly approximated elastic net
//'@description Object for smoothly approximated elastic net optimization with
//'bfgs optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a lambda and an alpha value.
//'@returns a list with fit results
typedef bfgsEnet<SEMCpp> bfgsEnetSEM;
RCPP_EXPOSED_CLASS_NODECL(bfgsEnetSEM)
  RCPP_MODULE(bfgsEnetSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<bfgsEnetSEM>( "bfgsEnetSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaEnet.")
    // methods
    .method( "setHessian", &bfgsEnetSEM::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &bfgsEnetSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

//'@name bfgsEnetMgSEM
//'@title smoothly approximated elastic net
//'@description Object for smoothly approximated elastic net optimization with
//'bfgs optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a lambda and an alpha value.
//'@returns a list with fit results
typedef bfgsEnet<mgSEM> bfgsEnetMgSEM;
RCPP_EXPOSED_CLASS_NODECL(bfgsEnetMgSEM)
  RCPP_MODULE(bfgsEnetMgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<bfgsEnetMgSEM>( "bfgsEnetMgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaEnet.")
    // methods
    .method( "setHessian", &bfgsEnetMgSEM::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &bfgsEnetMgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

//'@name glmnetEnetSEM
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
typedef glmnetEnet<SEMCpp> glmnetEnetSEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetEnetSEM)
  RCPP_MODULE(glmnetEnetSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetEnetSEM>( "glmnetEnetSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new glmnetEnetSEM.")
    // methods
    .method( "setHessian", &glmnetEnetSEM::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &glmnetEnetSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

//'@name glmnetEnetMgSEM
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
typedef glmnetEnet<mgSEM> glmnetEnetMgSEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetEnetMgSEM)
  RCPP_MODULE(glmnetEnetMgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetEnetMgSEM>( "glmnetEnetMgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new glmnetEnetMgSEM.")
    // methods
    .method( "setHessian", &glmnetEnetMgSEM::setHessian, "Changes the initial hessian. Expects a matrix")
    .method( "optimize", &glmnetEnetMgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

//'@name istaEnetSEM
//'@title elastic net optimization with ista optimizer
//'@description Object for elastic net optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a lambda and an alpha value.
//'@returns a list with fit results
//'
typedef istaEnet<SEMCpp> istaEnetSEM;
RCPP_EXPOSED_CLASS_NODECL(istaEnetSEM)
  RCPP_MODULE(istaEnetSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaEnetSEM>( "istaEnetSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaEnetSEM.")
    // methods
    .method( "optimize", &istaEnetSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }

//'@name istaEnetMgSEM
//'@title elastic net optimization with ista optimizer
//'@description Object for elastic net optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a lambda and an alpha value.
//'@returns a list with fit results
//'
typedef istaEnet<mgSEM> istaEnetMgSEM;
RCPP_EXPOSED_CLASS_NODECL(istaEnetMgSEM)
  RCPP_MODULE(istaEnetMgSEM_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaEnetMgSEM>( "istaEnetMgSEM" )
      .constructor<arma::rowvec,Rcpp::List>("Creates a new istaEnetMgSEM.")
    // methods
    .method( "optimize", &istaEnetMgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
    ;
  }
