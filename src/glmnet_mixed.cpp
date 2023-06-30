#include <RcppArmadillo.h>
#include "SEM.h"
#include "mgSEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

template<typename sem>
  class glmnetMixedPenalty{
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
    glmnetMixedPenalty(
      const arma::rowvec weights_,
      const std::vector<int> pType_int,
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
      
      pType.resize(pType_int.size());
      
      for(unsigned int i = 0; i < pType_int.size(); i++){
        pType.at(i) = static_cast<lessSEM::penaltyType>(pType_int.at(i));
      }
      
    }
    
    void setHessian(arma::mat newHessian){
      initialHessian = newHessian;
    }
    
    Rcpp::List optimize(Rcpp::NumericVector startingValues_, 
                        sem& SEM_, 
                        arma::rowvec lambda_, 
                        arma::rowvec theta_,
                        arma::rowvec alpha_){
      
      const double N = SEM_.sampleSize;
      
      SEMFitFramework<sem> SEMFF(SEM_, 1.0/N);
      
      lessSEM::tuningParametersMixedGlmnet tp;
      
      tp.penaltyType_ = pType;
      tp.lambda = lambda_;
      tp.theta = theta_;
      tp.alpha = alpha_;
      tp.weights = weights;
      
      lessSEM::penaltyMixedGlmnet penalty_;
      lessSEM::noSmoothPenalty<lessSEM::tuningParametersMixedGlmnet> smoothPenalty_;
      
      initializeMixedPenaltiesGlmnet(penalty_,
                                     pType);
      
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

//'@name glmnetMixedSEM
//'@title mixed optimization with glmnet optimizer
//'@description Object for mixed optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef glmnetMixedPenalty<SEMCpp> glmnetMixedSEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetMixedSEM)
RCPP_MODULE(glmnetMixedSEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<glmnetMixedSEM>( "glmnetMixedSEM" )
  .constructor<arma::rowvec,std::vector<int>,Rcpp::List>("Creates a new glmnetMixedSEM.")
  // methods
  .method( "setHessian", &glmnetMixedSEM::setHessian, "Changes the initial hessian. Expects a matrix")
  .method( "optimize", &glmnetMixedSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
  ;
}

//'@name glmnetMixedMgSEM
//'@title mixed optimization with glmnet optimizer
//'@description Object for mixed optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef glmnetMixedPenalty<mgSEM> glmnetMixedMgSEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetMixedMgSEM)
RCPP_MODULE(glmnetMixedMgSEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<glmnetMixedMgSEM>( "glmnetMixedMgSEM" )
  .constructor<arma::rowvec,std::vector<int>,Rcpp::List>("Creates a new glmnetMixedMgSEM.")
  // methods
  .method( "setHessian", &glmnetMixedMgSEM::setHessian, "Changes the initial hessian. Expects a matrix")
  .method( "optimize", &glmnetMixedMgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
  ;
}
