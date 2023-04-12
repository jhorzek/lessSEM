#include <RcppArmadillo.h>
#include "SEM.h"
#include "mgSEM.h"
#include "SEMFitFramework.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

template<typename sem>
  class glmnetCappedL1{
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
    glmnetCappedL1(
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
      
    }
    
    void setHessian(arma::mat newHessian){
      initialHessian = newHessian;
    }
    
    Rcpp::List optimize(Rcpp::NumericVector startingValues_, 
                        sem& SEM_, 
                        double theta_,
                        double lambda_){
      
      const double N = SEM_.sampleSize;
      
      SEMFitFramework<sem> SEMFF(SEM_);
      
      lessSEM::tuningParametersCappedL1Glmnet tp;
      
      tp.lambda = N*lambda_;
      tp.theta = theta_;
      tp.weights = weights;
      
      lessSEM::penaltyCappedL1Glmnet penalty_;
      lessSEM::noSmoothPenalty<lessSEM::tuningParametersCappedL1Glmnet> smoothPenalty_;
      
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

//'@name glmnetCappedL1SEM
//'@title CappedL1 optimization with glmnet optimizer
//'@description Object for cappedL1 optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef glmnetCappedL1<SEMCpp> glmnetCappedL1SEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetCappedL1SEM)
RCPP_MODULE(glmnetCappedL1SEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<glmnetCappedL1SEM>( "glmnetCappedL1SEM" )
  .constructor<arma::rowvec,Rcpp::List>("Creates a new glmnetCappedL1SEM.")
  // methods
  .method( "setHessian", &glmnetCappedL1SEM::setHessian, "Changes the initial hessian. Expects a matrix")
  .method( "optimize", &glmnetCappedL1SEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
  ;
}

//'@name glmnetCappedL1MgSEM
//'@title CappedL1 optimization with glmnet optimizer
//'@description Object for cappedL1 optimization with
//'glmnet optimizer
//'@field new creates a new object. Requires (2) a list with control elements
//'@field setHessian changes the Hessian of the model. Expects a matrix
//'@field optimize optimize the model. Expects a vector with starting values,
//'a SEM of type SEM_Cpp, a theta and a lambda value.
//'@returns a list with fit results
//'
typedef glmnetCappedL1<mgSEM> glmnetCappedL1MgSEM;
RCPP_EXPOSED_CLASS_NODECL(glmnetCappedL1MgSEM)
RCPP_MODULE(glmnetCappedL1MgSEM_cpp){
  using namespace Rcpp;
  Rcpp::class_<glmnetCappedL1MgSEM>( "glmnetCappedL1MgSEM" )
  .constructor<arma::rowvec,Rcpp::List>("Creates a new glmnetCappedL1MgSEM.")
  // methods
  .method( "setHessian", &glmnetCappedL1MgSEM::setHessian, "Changes the initial hessian. Expects a matrix")
  .method( "optimize", &glmnetCappedL1MgSEM::optimize, "Optimizes the model. Expects SEM, labeled vector with starting values, lambda, and alpha")
  ;
}
