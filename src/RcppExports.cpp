// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// innerGLMNET
arma::colvec innerGLMNET(const arma::colvec parameters, const int N, const arma::colvec subGroupGradient, const arma::mat subGroupHessian, const double subGroupLambda, const Rcpp::LogicalVector regularized, const arma::colvec adaptiveLassoWeights, const int maxIter, const double epsBreak, const bool useMultipleConvergenceCriteria);
RcppExport SEXP _linr_innerGLMNET(SEXP parametersSEXP, SEXP NSEXP, SEXP subGroupGradientSEXP, SEXP subGroupHessianSEXP, SEXP subGroupLambdaSEXP, SEXP regularizedSEXP, SEXP adaptiveLassoWeightsSEXP, SEXP maxIterSEXP, SEXP epsBreakSEXP, SEXP useMultipleConvergenceCriteriaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type subGroupGradient(subGroupGradientSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type subGroupHessian(subGroupHessianSEXP);
    Rcpp::traits::input_parameter< const double >::type subGroupLambda(subGroupLambdaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type regularized(regularizedSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type adaptiveLassoWeights(adaptiveLassoWeightsSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type epsBreak(epsBreakSEXP);
    Rcpp::traits::input_parameter< const bool >::type useMultipleConvergenceCriteria(useMultipleConvergenceCriteriaSEXP);
    rcpp_result_gen = Rcpp::wrap(innerGLMNET(parameters, N, subGroupGradient, subGroupHessian, subGroupLambda, regularized, adaptiveLassoWeights, maxIter, epsBreak, useMultipleConvergenceCriteria));
    return rcpp_result_gen;
END_RCPP
}
// computeIndividualM2LL
double computeIndividualM2LL(const int nObservedVariables, const arma::colvec& rawData, const arma::colvec& impliedMeans, const arma::mat& impliedCovariance);
RcppExport SEXP _linr_computeIndividualM2LL(SEXP nObservedVariablesSEXP, SEXP rawDataSEXP, SEXP impliedMeansSEXP, SEXP impliedCovarianceSEXP) {
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
RcppExport SEXP _linr_computeGroupM2LL(SEXP sampleSizeSEXP, SEXP nObservedVariablesSEXP, SEXP observedMeansSEXP, SEXP observedCovSEXP, SEXP impliedMeansSEXP, SEXP impliedCovarianceSEXP) {
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
RcppExport SEXP _linr_computeImpliedCovariance(SEXP FmatrixSEXP, SEXP AmatrixSEXP, SEXP SmatrixSEXP) {
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
arma::mat computeImpliedMeans(const arma::mat& Fmatrix, const arma::mat& Amatrix, const arma::colvec& Mvector);
RcppExport SEXP _linr_computeImpliedMeans(SEXP FmatrixSEXP, SEXP AmatrixSEXP, SEXP MvectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Fmatrix(FmatrixSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Amatrix(AmatrixSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Mvector(MvectorSEXP);
    rcpp_result_gen = Rcpp::wrap(computeImpliedMeans(Fmatrix, Amatrix, Mvector));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_SEM_cpp();
RcppExport SEXP _rcpp_module_boot_glmnetEnet_cpp();
RcppExport SEXP _rcpp_module_boot_istaEnet_cpp();
RcppExport SEXP _rcpp_module_boot_istaEnetGeneralPurpose_cpp();

static const R_CallMethodDef CallEntries[] = {
    {"_linr_innerGLMNET", (DL_FUNC) &_linr_innerGLMNET, 10},
    {"_linr_computeIndividualM2LL", (DL_FUNC) &_linr_computeIndividualM2LL, 4},
    {"_linr_computeGroupM2LL", (DL_FUNC) &_linr_computeGroupM2LL, 6},
    {"_linr_computeImpliedCovariance", (DL_FUNC) &_linr_computeImpliedCovariance, 3},
    {"_linr_computeImpliedMeans", (DL_FUNC) &_linr_computeImpliedMeans, 3},
    {"_rcpp_module_boot_SEM_cpp", (DL_FUNC) &_rcpp_module_boot_SEM_cpp, 0},
    {"_rcpp_module_boot_glmnetEnet_cpp", (DL_FUNC) &_rcpp_module_boot_glmnetEnet_cpp, 0},
    {"_rcpp_module_boot_istaEnet_cpp", (DL_FUNC) &_rcpp_module_boot_istaEnet_cpp, 0},
    {"_rcpp_module_boot_istaEnetGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaEnetGeneralPurpose_cpp, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_linr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
