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
RcppExport SEXP _lessSEM_innerGLMNET(SEXP parametersSEXP, SEXP NSEXP, SEXP subGroupGradientSEXP, SEXP subGroupHessianSEXP, SEXP subGroupLambdaSEXP, SEXP regularizedSEXP, SEXP adaptiveLassoWeightsSEXP, SEXP maxIterSEXP, SEXP epsBreakSEXP, SEXP useMultipleConvergenceCriteriaSEXP) {
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
// glmnetInner_C
arma::rowvec glmnetInner_C(const arma::rowvec& parameters_kMinus1, const arma::rowvec& gradients_kMinus1, const arma::mat& Hessian, const double lambda, const double alpha, const arma::rowvec& weights, const int maxIterIn, const double breakInner, const int verbose);
RcppExport SEXP _lessSEM_glmnetInner_C(SEXP parameters_kMinus1SEXP, SEXP gradients_kMinus1SEXP, SEXP HessianSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP weightsSEXP, SEXP maxIterInSEXP, SEXP breakInnerSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type parameters_kMinus1(parameters_kMinus1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type gradients_kMinus1(gradients_kMinus1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hessian(HessianSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIterIn(maxIterInSEXP);
    Rcpp::traits::input_parameter< const double >::type breakInner(breakInnerSEXP);
    Rcpp::traits::input_parameter< const int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(glmnetInner_C(parameters_kMinus1, gradients_kMinus1, Hessian, lambda, alpha, weights, maxIterIn, breakInner, verbose));
    return rcpp_result_gen;
END_RCPP
}
// BFGS_C
arma::mat BFGS_C(const arma::rowvec& parameters_kMinus1, const arma::rowvec& gradients_kMinus1, const arma::mat& Hessian_kMinus1, const arma::rowvec& parameters_k, const arma::rowvec& gradients_k, const bool cautious, const double hessianEps);
RcppExport SEXP _lessSEM_BFGS_C(SEXP parameters_kMinus1SEXP, SEXP gradients_kMinus1SEXP, SEXP Hessian_kMinus1SEXP, SEXP parameters_kSEXP, SEXP gradients_kSEXP, SEXP cautiousSEXP, SEXP hessianEpsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type parameters_kMinus1(parameters_kMinus1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type gradients_kMinus1(gradients_kMinus1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Hessian_kMinus1(Hessian_kMinus1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type parameters_k(parameters_kSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type gradients_k(gradients_kSEXP);
    Rcpp::traits::input_parameter< const bool >::type cautious(cautiousSEXP);
    Rcpp::traits::input_parameter< const double >::type hessianEps(hessianEpsSEXP);
    rcpp_result_gen = Rcpp::wrap(BFGS_C(parameters_kMinus1, gradients_kMinus1, Hessian_kMinus1, parameters_k, gradients_k, cautious, hessianEps));
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
RcppExport SEXP _rcpp_module_boot_glmnetEnetGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaEnetGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaLspGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaMcpGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaScadGeneralPurpose_cpp();
RcppExport SEXP _rcpp_module_boot_istaLSP_cpp();
RcppExport SEXP _rcpp_module_boot_istaMcp_cpp();
RcppExport SEXP _rcpp_module_boot_istaScad_cpp();

static const R_CallMethodDef CallEntries[] = {
    {"_lessSEM_innerGLMNET", (DL_FUNC) &_lessSEM_innerGLMNET, 10},
    {"_lessSEM_glmnetInner_C", (DL_FUNC) &_lessSEM_glmnetInner_C, 9},
    {"_lessSEM_BFGS_C", (DL_FUNC) &_lessSEM_BFGS_C, 7},
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
    {"_rcpp_module_boot_glmnetEnetGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_glmnetEnetGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaEnetGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaEnetGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaLspGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaLspGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaMcpGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaMcpGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaScadGeneralPurpose_cpp", (DL_FUNC) &_rcpp_module_boot_istaScadGeneralPurpose_cpp, 0},
    {"_rcpp_module_boot_istaLSP_cpp", (DL_FUNC) &_rcpp_module_boot_istaLSP_cpp, 0},
    {"_rcpp_module_boot_istaMcp_cpp", (DL_FUNC) &_rcpp_module_boot_istaMcp_cpp, 0},
    {"_rcpp_module_boot_istaScad_cpp", (DL_FUNC) &_rcpp_module_boot_istaScad_cpp, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_lessSEM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
