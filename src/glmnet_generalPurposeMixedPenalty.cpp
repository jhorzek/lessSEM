#include <RcppArmadillo.h>
#include "SEM.h"
#include "gpFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

  class glmnetMixedPenaltyGeneralPurpose{
    public:
      
    std::vector<lessSEM::penaltyType> pType;
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
    glmnetMixedPenaltyGeneralPurpose(
      const arma::rowvec weights_,
      const std::vector<std::string> pType_str,
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
      
      pType = lessSEM::stringPenaltyToPenaltyType(pType_str);
      
    }
    
    void setHessian(arma::mat newHessian){
      initialHessian = newHessian;
    }
    
    Rcpp::List optimize(Rcpp::NumericVector startingValues_, 
                        Rcpp::Function fitFunction,
                        Rcpp::Function gradientFunction,
                        Rcpp::List userSuppliedElements,
                        arma::rowvec lambda_, 
                        arma::rowvec theta_,
                        arma::rowvec alpha_){
      
      generalPurposeFitFramework gpFF(fitFunction, 
                                      gradientFunction, 
                                      userSuppliedElements);
      
      const double N = 1;
      
      lessSEM::tuningParametersMixedGlmnet tp;
      
      tp.penaltyType_ = pType;
      tp.lambda = lambda_;
      tp.theta = theta_;
      tp.alpha = alpha_;
      tp.weights = weights;
      
      lessSEM::penaltyMixedGlmnet penalty_;
      lessSEM::noSmoothPenalty<lessSEM::tuningParametersMixedGlmnet> smoothPenalty_;
      
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
        gpFF,
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

  class glmnetMixedPenaltyGeneralPurposeCpp{
    public:
      
      std::vector<lessSEM::penaltyType> pType;
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
    glmnetMixedPenaltyGeneralPurposeCpp(
      const arma::rowvec weights_,
      const std::vector<std::string> pType_str,
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
      
      pType = lessSEM::stringPenaltyToPenaltyType(pType_str);
      
    }
    
    void setHessian(arma::mat newHessian){
      initialHessian = newHessian;
    }
    
    Rcpp::List optimize(Rcpp::NumericVector startingValues_, 
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
      
      const double N = 1;
      
      lessSEM::tuningParametersMixedGlmnet tp;
      
      tp.penaltyType_ = pType;
      tp.lambda = lambda_;
      tp.theta = theta_;
      tp.alpha = alpha_;
      tp.weights = weights;
      
      lessSEM::penaltyMixedGlmnet penalty_;
      lessSEM::noSmoothPenalty<lessSEM::tuningParametersMixedGlmnet> smoothPenalty_;
      
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
        gpFF,
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

//'@name glmnetMixedPenaltyGeneralPurpose
//'@title mixed optimization with glmnet optimizer
//'@description Object for mixed optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
RCPP_EXPOSED_CLASS(glmnetMixedPenaltyGeneralPurpose)
RCPP_MODULE(glmnetMixedPenaltyGeneralPurpose_cpp){
  using namespace Rcpp;
  Rcpp::class_<glmnetMixedPenaltyGeneralPurpose>( "glmnetMixedPenaltyGeneralPurpose" )
  .constructor<arma::rowvec, std::vector<std::string>,Rcpp::List>("Creates a new glmnetMixedPenaltyGeneralPurpose")
  // methods
  .method( "setHessian", &glmnetMixedPenaltyGeneralPurpose::setHessian, "Changes the initial hessian. Expects a matrix")
  .method( "optimize", &glmnetMixedPenaltyGeneralPurpose::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values, theta, and lambda")
  ;
}

//'@name glmnetMixedPenaltyGeneralPurposeCpp
//'@title mixed optimization with glmnet optimizer
//'@description Object for mixed optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
RCPP_EXPOSED_CLASS(glmnetMixedPenaltyGeneralPurposeCpp)
RCPP_MODULE(glmnetMixedPenaltyGeneralPurposeCpp_cpp){
  using namespace Rcpp;
  Rcpp::class_<glmnetMixedPenaltyGeneralPurposeCpp>( "glmnetMixedPenaltyGeneralPurposeCpp" )
  .constructor<arma::rowvec, std::vector<std::string>,Rcpp::List>("Creates a new glmnetMixedPenaltyGeneralPurposeCpp")
  // methods
  .method( "setHessian", &glmnetMixedPenaltyGeneralPurposeCpp::setHessian, "Changes the initial hessian. Expects a matrix")
  .method( "optimize", &glmnetMixedPenaltyGeneralPurposeCpp::optimize, "Optimizes the model. Expects fitFunction, gradientFunction, userSuppliedElements, labeled vector with starting values, theta, and lambda")
  ;
}