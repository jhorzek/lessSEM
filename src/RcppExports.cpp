// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/lessSEM.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// callFitFunction
double callFitFunction(SEXP fitFunctionSEXP, Rcpp::NumericVector parameters, Rcpp::List userSuppliedElements);
RcppExport SEXP _lessSEM_callFitFunction(SEXP fitFunctionSEXPSEXP, SEXP parametersSEXP, SEXP userSuppliedElementsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type fitFunctionSEXP(fitFunctionSEXPSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type userSuppliedElements(userSuppliedElementsSEXP);
    rcpp_result_gen = Rcpp::wrap(callFitFunction(fitFunctionSEXP, parameters, userSuppliedElements));
    return rcpp_result_gen;
END_RCPP
}
// computeIndividualM2LL
double computeIndividualM2LL(const int nObservedVariables, const arma::colvec& rawData, const arma::colvec& impliedMeans, const arma::mat& impliedCovariance);
RcppExport SEXP _lessSEM_computeIndividualM2LL(SEXP nObservedVariablesSEXP, SEXP rawDataSEXP, SEXP impliedMeansSEXP, SEXP impliedCovarianceSEXP) {
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
RcppExport SEXP _lessSEM_computeGroupM2LL(SEXP sampleSizeSEXP, SEXP nObservedVariablesSEXP, SEXP observedMeansSEXP, SEXP observedCovSEXP, SEXP impliedMeansSEXP, SEXP impliedCovarianceSEXP) {
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
RcppExport SEXP _lessSEM_computeImpliedCovariance(SEXP FmatrixSEXP, SEXP AmatrixSEXP, SEXP SmatrixSEXP) {
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
RcppExport SEXP _lessSEM_computeImpliedMeans(SEXP FmatrixSEXP, SEXP AmatrixSEXP, SEXP MvectorSEXP) {
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
// mcpPenalty_C
double mcpPenalty_C(const double par, const double lambda_p, const double theta);
RcppExport SEXP _lessSEM_mcpPenalty_C(SEXP parSEXP, SEXP lambda_pSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type par(parSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_p(lambda_pSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(mcpPenalty_C(par, lambda_p, theta));
    return rcpp_result_gen;
END_RCPP
}
// scadPenalty_C
double scadPenalty_C(const double par, const double lambda_p, const double theta);
RcppExport SEXP _lessSEM_scadPenalty_C(SEXP parSEXP, SEXP lambda_pSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type par(parSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_p(lambda_pSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(scadPenalty_C(par, lambda_p, theta));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_SEM_cpp();
RcppExport SEXP _rcpp_module_boot_istaCappedL1_cpp();
RcppExport SEXP _rcpp_module_boot_bfgsEnet_cpp();
RcppExport SEXP _rcpp_module_boot_glmnetEnet_cpp();
RcppExport SEXP _rcpp_module_boot_istaEnet_cpp();
RcppExport SEXP _rcpp_module_boot_istaCappedL1GeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaCappedL1GeneralPurposeCpp_cpp();
RcppExport SEXP _rcpp_module_boot_glmnetEnetGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaEnetGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_glmnetEnetGeneralPurposeCpp_cpp();
RcppExport SEXP _rcpp_module_boot_istaEnetGeneralPurposeCpp_cpp();
RcppExport SEXP _rcpp_module_boot_istaLspGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaLspGeneralPurposeCpp_cpp();
RcppExport SEXP _rcpp_module_boot_istaMcpGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaMcpGeneralPurposeCpp_cpp();
RcppExport SEXP _rcpp_module_boot_istaScadGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaScadGeneralPurposeCpp_cpp();
RcppExport SEXP _rcpp_module_boot_istaLSP_cpp();
RcppExport SEXP _rcpp_module_boot_istaMcp_cpp();
RcppExport SEXP _rcpp_module_boot_mgSEM_cpp();
RcppExport SEXP _rcpp_module_boot_istaScad_cpp();

static const R_CallMethodDef CallEntries[] = {
    {"_lessSEM_callFitFunction", (DL_FUNC) &_lessSEM_callFitFunction, 3},
    {"_lessSEM_computeIndividualM2LL", (DL_FUNC) &_lessSEM_computeIndividualM2LL, 4},
    {"_lessSEM_computeGroupM2LL", (DL_FUNC) &_lessSEM_computeGroupM2LL, 6},
    {"_lessSEM_computeImpliedCovariance", (DL_FUNC) &_lessSEM_computeImpliedCovariance, 3},
    {"_lessSEM_computeImpliedMeans", (DL_FUNC) &_lessSEM_computeImpliedMeans, 3},
    {"_lessSEM_mcpPenalty_C", (DL_FUNC) &_lessSEM_mcpPenalty_C, 3},
    {"_lessSEM_scadPenalty_C", (DL_FUNC) &_lessSEM_scadPenalty_C, 3},
    {"_rcpp_module_boot_SEM_cpp", (DL_FUNC) &_rcpp_module_boot_SEM_cpp, 0},
    {"_rcpp_module_boot_istaCappedL1_cpp", (DL_FUNC) &_rcpp_module_boot_istaCappedL1_cpp, 0},
    {"_rcpp_module_boot_bfgsEnet_cpp", (DL_FUNC) &_rcpp_module_boot_bfgsEnet_cpp, 0},
    {"_rcpp_module_boot_glmnetEnet_cpp", (DL_FUNC) &_rcpp_module_boot_glmnetEnet_cpp, 0},
    {"_rcpp_module_boot_istaEnet_cpp", (DL_FUNC) &_rcpp_module_boot_istaEnet_cpp, 0},
    {"_rcpp_module_boot_istaCappedL1GeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaCappedL1GeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaCappedL1GeneralPurposeCpp_cpp", (DL_FUNC) &_rcpp_module_boot_istaCappedL1GeneralPurposeCpp_cpp, 0},
    {"_rcpp_module_boot_glmnetEnetGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_glmnetEnetGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaEnetGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaEnetGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_glmnetEnetGeneralPurposeCpp_cpp", (DL_FUNC) &_rcpp_module_boot_glmnetEnetGeneralPurposeCpp_cpp, 0},
    {"_rcpp_module_boot_istaEnetGeneralPurposeCpp_cpp", (DL_FUNC) &_rcpp_module_boot_istaEnetGeneralPurposeCpp_cpp, 0},
    {"_rcpp_module_boot_istaLspGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaLspGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaLspGeneralPurposeCpp_cpp", (DL_FUNC) &_rcpp_module_boot_istaLspGeneralPurposeCpp_cpp, 0},
    {"_rcpp_module_boot_istaMcpGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaMcpGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaMcpGeneralPurposeCpp_cpp", (DL_FUNC) &_rcpp_module_boot_istaMcpGeneralPurposeCpp_cpp, 0},
    {"_rcpp_module_boot_istaScadGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaScadGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaScadGeneralPurposeCpp_cpp", (DL_FUNC) &_rcpp_module_boot_istaScadGeneralPurposeCpp_cpp, 0},
    {"_rcpp_module_boot_istaLSP_cpp", (DL_FUNC) &_rcpp_module_boot_istaLSP_cpp, 0},
    {"_rcpp_module_boot_istaMcp_cpp", (DL_FUNC) &_rcpp_module_boot_istaMcp_cpp, 0},
    {"_rcpp_module_boot_mgSEM_cpp", (DL_FUNC) &_rcpp_module_boot_mgSEM_cpp, 0},
    {"_rcpp_module_boot_istaScad_cpp", (DL_FUNC) &_rcpp_module_boot_istaScad_cpp, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_lessSEM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
