#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

// creating types for user supplied C++ functions. See
// https://gallery.rcpp.org/articles/passing-cpp-function-pointers/
typedef double (*fitFunPtr)(const Rcpp::NumericVector&, //parameters
                Rcpp::List&); //additional elements
typedef Rcpp::XPtr<fitFunPtr> fitFunPtr_t;

//' callFitFunction
//' 
//' wrapper to call user defined fit function
//' @param fitFunctionSEXP pointer to fit function
//' @param parameters vector with parameter values
//' @param userSuppliedElements list with additional elements
//' @returns fit value (double)
//' 
// [[Rcpp::export]]
double callFitFunction(SEXP fitFunctionSEXP, 
                       Rcpp::NumericVector parameters, 
                       Rcpp::List userSuppliedElements){
  
  Rcpp::XPtr<fitFunPtr> xpFitFun(fitFunctionSEXP);
  fitFunPtr fitFunction = *xpFitFun;
  
  return(fitFunction(parameters,
                     userSuppliedElements));
  
}






