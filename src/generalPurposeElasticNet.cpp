#include <RcppArmadillo.h>
#include "SEM.h"
#include "gpFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

//'@name istaEnetGeneralPurpose
//'@title elastic net optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'an R function to compute the fit, an R function to compute the gradients, a
//'list with elements the fit and gradient function require, a lambda and an alpha value.
//'@returns a list with fit results
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
  const lessSEM::convCritInnerIsta convCritInner;
  const double sigma;
  const lessSEM::stepSizeInheritance stepSizeInh;
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
    convCritInner(static_cast<lessSEM::convCritInnerIsta>(Rcpp::as<int> (control["convCritInner"]))),
    sigma(Rcpp::as<double> (control["sigma"])),
    stepSizeInh(static_cast<lessSEM::stepSizeInheritance>(Rcpp::as<int> (control["stepSizeInheritance"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      Rcpp::Function fitFunction,
      Rcpp::Function gradientFunction,
      Rcpp::List userSuppliedElements,
      double lambda_,
      double alpha_){
    
    generalPurposeFitFramework gpFF(fitFunction, gradientFunction, userSuppliedElements);
    
    lessSEM::tuningParametersEnet tp;
    tp.lambda = lambda_;
    tp.weights = weights;
    tp.alpha = alpha_;
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
      1,
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::ista(
      gpFF,
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

//'@name glmnetEnetGeneralPurpose
//'@title elastic net optimization with glmnet optimizer
//'@description Object for elastic net optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'an R function to compute the fit, an R function to compute the gradients, a
//'list with elements the fit and gradient function require, a lambda and an alpha value.
//'@returns a list with fit results
//'
class glmnetEnetGeneralPurpose{
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
    
    convergenceCriterion(static_cast<lessSEM::convergenceCriteriaGlmnet>(Rcpp::as<int> (control["convergenceCriterion"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  void setHessian(arma::mat newHessian){
    initialHessian = newHessian;
  }
  
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
    
    lessSEM::tuningParametersEnet tp;
    tp.lambda = lambda_;
    tp.weights = weights;
    tp.alpha = alpha_;
    
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
      breakOuter, // change in fit required to break the outer iteration
      breakInner,
      convergenceCriterion, // this is related to the inner
      // breaking condition. 
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::glmnet(
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

//'@name istaEnetGeneralPurposeCpp
//'@title elastic net optimization with ista
//'@description Object for elastic net optimization with
//'ista optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEXP function pointer to compute the fit, a SEXP function pointer to compute the gradients, a
//'list with elements the fit and gradient function require, a lambda and an alpha value.
//'@returns a list with fit results
class istaEnetGeneralPurposeCpp{
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
  const lessSEM::convCritInnerIsta convCritInner;
  const double sigma;
  const lessSEM::stepSizeInheritance stepSizeInh;
  const int verbose; 
  
  // constructor
  istaEnetGeneralPurposeCpp(
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
    convCritInner(static_cast<lessSEM::convCritInnerIsta>(Rcpp::as<int> (control["convCritInner"]))),
    sigma(Rcpp::as<double> (control["sigma"])),
    stepSizeInh(static_cast<lessSEM::stepSizeInheritance>(Rcpp::as<int> (control["stepSizeInheritance"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      SEXP fitFunction,
      SEXP gradientFunction,
      Rcpp::List userSuppliedElements,
      double lambda_,
      double alpha_){
    
    generalPurposeFitFrameworkCpp gpFF(startingValues_, fitFunction, gradientFunction, userSuppliedElements);
    
    lessSEM::tuningParametersEnet tp;
    tp.lambda = lambda_;
    tp.weights = weights;
    tp.alpha = alpha_;
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
      1,
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::ista(
      gpFF,
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

//'@name glmnetEnetGeneralPurposeCpp
//'@title elastic net optimization with glmnet optimizer
//'@description Object for elastic net optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (1) a vector with weights for each
//'parameter and (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEXP function pointer to compute the fit, a SEXP function pointer to compute the gradients, a
//'list with elements the fit and gradient function require, a lambda and an alpha value.
//'@returns a list with fit results
class glmnetEnetGeneralPurposeCpp{
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
  glmnetEnetGeneralPurposeCpp(
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
    
    convergenceCriterion(static_cast<lessSEM::convergenceCriteriaGlmnet>(Rcpp::as<int> (control["convergenceCriterion"]))),
    verbose(Rcpp::as<int> (control["verbose"])){}
  
  void setHessian(arma::mat newHessian){
    initialHessian = newHessian;
  }
  
  Rcpp::List optimize(
      Rcpp::NumericVector startingValues_, 
      SEXP fitFunction,
      SEXP gradientFunction,
      Rcpp::List userSuppliedElements,
      double lambda_,
      double alpha_){
    
    generalPurposeFitFrameworkCpp gpFF(startingValues_, fitFunction, gradientFunction, userSuppliedElements);
    
    lessSEM::tuningParametersEnet tp;
    tp.lambda = lambda_;
    tp.weights = weights;
    tp.alpha = alpha_;
    
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
      breakOuter, // change in fit required to break the outer iteration
      breakInner,
      convergenceCriterion, // this is related to the inner
      // breaking condition. 
      verbose
    };
    
    lessSEM::fitResults fitResults_ = lessSEM::glmnet(
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

/////////////

RCPP_EXPOSED_CLASS(glmnetEnetGeneralPurpose)
  RCPP_MODULE(glmnetEnetGeneralPurpose_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetEnetGeneralPurpose>( "glmnetEnetGeneralPurpose" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new glmnetEnetGeneralPurpose.")
    // methods
    .method( "optimize", &glmnetEnetGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values and lambda")
    .method( "setHessian", &glmnetEnetGeneralPurpose::setHessian, "Change the initial Hessian matrix.")
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

RCPP_EXPOSED_CLASS(glmnetEnetGeneralPurposeCpp)
  RCPP_MODULE(glmnetEnetGeneralPurposeCpp_cpp){
    using namespace Rcpp;
    Rcpp::class_<glmnetEnetGeneralPurposeCpp>( "glmnetEnetGeneralPurposeCpp" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new glmnetEnetGeneralPurposeCpp.")
    // methods
    .method( "optimize", &glmnetEnetGeneralPurposeCpp::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values and lambda")
    .method( "setHessian", &glmnetEnetGeneralPurposeCpp::setHessian, "Change the initial Hessian matrix.")
    ;
  }

RCPP_EXPOSED_CLASS(istaEnetGeneralPurposeCpp)
  RCPP_MODULE(istaEnetGeneralPurposeCpp_cpp){
    using namespace Rcpp;
    Rcpp::class_<istaEnetGeneralPurposeCpp>( "istaEnetGeneralPurposeCpp" )
      .constructor<Rcpp::NumericVector,Rcpp::List>("Creates a new istaEnetGeneralPurposeCpp.")
    // methods
    .method( "optimize", &istaEnetGeneralPurposeCpp::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values and lambda")
    ;
  }