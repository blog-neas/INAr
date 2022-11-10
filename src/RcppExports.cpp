// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sunMC_Cpp
NumericVector sunMC_Cpp(NumericVector x, unsigned int method);
RcppExport SEXP _INAr_sunMC_Cpp(SEXP xSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(sunMC_Cpp(x, method));
    return rcpp_result_gen;
END_RCPP
}
// sunMC_semiparBOOT_Cpp
NumericVector sunMC_semiparBOOT_Cpp(NumericVector x, int B, unsigned int method);
RcppExport SEXP _INAr_sunMC_semiparBOOT_Cpp(SEXP xSEXP, SEXP BSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(sunMC_semiparBOOT_Cpp(x, B, method));
    return rcpp_result_gen;
END_RCPP
}
// sunMC_parBOOT_Cpp
NumericVector sunMC_parBOOT_Cpp(NumericVector x, int B, unsigned int method);
RcppExport SEXP _INAr_sunMC_parBOOT_Cpp(SEXP xSEXP, SEXP BSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(sunMC_parBOOT_Cpp(x, B, method));
    return rcpp_result_gen;
END_RCPP
}
// sunMC_pitBOOT_Cpp
NumericVector sunMC_pitBOOT_Cpp(NumericVector x, int B, unsigned int method);
RcppExport SEXP _INAr_sunMC_pitBOOT_Cpp(SEXP xSEXP, SEXP BSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(sunMC_pitBOOT_Cpp(x, B, method));
    return rcpp_result_gen;
END_RCPP
}
// sunMCtest_boot
List sunMCtest_boot(NumericVector X, unsigned int arrival, unsigned int method, int B);
RcppExport SEXP _INAr_sunMCtest_boot(SEXP XSEXP, SEXP arrivalSEXP, SEXP methodSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type arrival(arrivalSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(sunMCtest_boot(X, arrival, method, B));
    return rcpp_result_gen;
END_RCPP
}
// INARp_cpp
NumericVector INARp_cpp(NumericVector resid, DoubleVector a);
RcppExport SEXP _INAr_INARp_cpp(SEXP residSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type resid(residSEXP);
    Rcpp::traits::input_parameter< DoubleVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(INARp_cpp(resid, a));
    return rcpp_result_gen;
END_RCPP
}
// Xresid
Rcpp::List Xresid(NumericVector X, NumericVector alphas, double mINN, double vINN);
RcppExport SEXP _INAr_Xresid(SEXP XSEXP, SEXP alphasSEXP, SEXP mINNSEXP, SEXP vINNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< double >::type mINN(mINNSEXP);
    Rcpp::traits::input_parameter< double >::type vINN(vINNSEXP);
    rcpp_result_gen = Rcpp::wrap(Xresid(X, alphas, mINN, vINN));
    return rcpp_result_gen;
END_RCPP
}
// Xmoments
Rcpp::List Xmoments(NumericVector X, NumericVector alphas);
RcppExport SEXP _INAr_Xmoments(SEXP XSEXP, SEXP alphasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphas(alphasSEXP);
    rcpp_result_gen = Rcpp::wrap(Xmoments(X, alphas));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_INAr_sunMC_Cpp", (DL_FUNC) &_INAr_sunMC_Cpp, 2},
    {"_INAr_sunMC_semiparBOOT_Cpp", (DL_FUNC) &_INAr_sunMC_semiparBOOT_Cpp, 3},
    {"_INAr_sunMC_parBOOT_Cpp", (DL_FUNC) &_INAr_sunMC_parBOOT_Cpp, 3},
    {"_INAr_sunMC_pitBOOT_Cpp", (DL_FUNC) &_INAr_sunMC_pitBOOT_Cpp, 3},
    {"_INAr_sunMCtest_boot", (DL_FUNC) &_INAr_sunMCtest_boot, 4},
    {"_INAr_INARp_cpp", (DL_FUNC) &_INAr_INARp_cpp, 2},
    {"_INAr_Xresid", (DL_FUNC) &_INAr_Xresid, 4},
    {"_INAr_Xmoments", (DL_FUNC) &_INAr_Xmoments, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_INAr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
