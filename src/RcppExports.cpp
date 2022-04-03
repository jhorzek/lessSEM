// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// computeIndividualM2LL
double computeIndividualM2LL(const int nObservedVariables, const arma::colvec& rawData, const arma::colvec& impliedMeans, const arma::mat& impliedCovariance);
RcppExport SEXP _aCV4SEM_computeIndividualM2LL(SEXP nObservedVariablesSEXP, SEXP rawDataSEXP, SEXP impliedMeansSEXP, SEXP impliedCovarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nObservedVariables(nObservedVariablesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type rawData(rawDataSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type impliedMeans(impliedMeansSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type impliedCovariance(impliedCovarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(computeIndividualM2LL(nObservedVariables, rawData, impliedMeans, impliedCovariance));
    return rcpp_result_gen;
END_RCPP
}
// computeGroupM2LL
double computeGroupM2LL(const int sampleSize, const int nObservedVariables, const arma::colvec& observedMeans, const arma::mat& observedCov, const arma::colvec& impliedMeans, const arma::mat& impliedCovariance);
RcppExport SEXP _aCV4SEM_computeGroupM2LL(SEXP sampleSizeSEXP, SEXP nObservedVariablesSEXP, SEXP observedMeansSEXP, SEXP observedCovSEXP, SEXP impliedMeansSEXP, SEXP impliedCovarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type sampleSize(sampleSizeSEXP);
    Rcpp::traits::input_parameter< const int >::type nObservedVariables(nObservedVariablesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type observedMeans(observedMeansSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type observedCov(observedCovSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type impliedMeans(impliedMeansSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type impliedCovariance(impliedCovarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(computeGroupM2LL(sampleSize, nObservedVariables, observedMeans, observedCov, impliedMeans, impliedCovariance));
    return rcpp_result_gen;
END_RCPP
}
// computeImpliedCovariance
arma::mat computeImpliedCovariance(const arma::mat& Fmatrix, const arma::mat& Amatrix, const arma::mat& Smatrix);
RcppExport SEXP _aCV4SEM_computeImpliedCovariance(SEXP FmatrixSEXP, SEXP AmatrixSEXP, SEXP SmatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Fmatrix(FmatrixSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Amatrix(AmatrixSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Smatrix(SmatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(computeImpliedCovariance(Fmatrix, Amatrix, Smatrix));
    return rcpp_result_gen;
END_RCPP
}
// computeImpliedMeans
arma::mat computeImpliedMeans(const arma::mat& Fmatrix, const arma::mat& Amatrix, const arma::colvec& mVector);
RcppExport SEXP _aCV4SEM_computeImpliedMeans(SEXP FmatrixSEXP, SEXP AmatrixSEXP, SEXP mVectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Fmatrix(FmatrixSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Amatrix(AmatrixSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type mVector(mVectorSEXP);
    rcpp_result_gen = Rcpp::wrap(computeImpliedMeans(Fmatrix, Amatrix, mVector));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_SEM_cpp();

static const R_CallMethodDef CallEntries[] = {
    {"_aCV4SEM_computeIndividualM2LL", (DL_FUNC) &_aCV4SEM_computeIndividualM2LL, 4},
    {"_aCV4SEM_computeGroupM2LL", (DL_FUNC) &_aCV4SEM_computeGroupM2LL, 6},
    {"_aCV4SEM_computeImpliedCovariance", (DL_FUNC) &_aCV4SEM_computeImpliedCovariance, 3},
    {"_aCV4SEM_computeImpliedMeans", (DL_FUNC) &_aCV4SEM_computeImpliedMeans, 3},
    {"_rcpp_module_boot_SEM_cpp", (DL_FUNC) &_rcpp_module_boot_SEM_cpp, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_aCV4SEM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
