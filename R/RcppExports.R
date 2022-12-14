# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Sun McCabe score statistic to test for dependence in an integer autoregressive process
#' @param x NumericVector
#' @param method unsigned int
#' @details
#' This is an internal function, it will be excluded in future versions.
#' @export
sunMC_Cpp <- function(x, method) {
    .Call('_INAr_sunMC_Cpp', PACKAGE = 'INAr', x, method)
}

#' Semiparametric bootstrap version of the Sun McCabe score test.
#' @param x NumericVector
#' @param B int
#' @param method unsigned int
#' @details
#' This is an internal function, it will be excluded in future versions.
#' @export
sunMC_semiparBOOT_Cpp <- function(x, B, method) {
    .Call('_INAr_sunMC_semiparBOOT_Cpp', PACKAGE = 'INAr', x, B, method)
}

#' Parametric bootstrap version of the Sun McCabe score test.
#' @param x NumericVector
#' @param B int
#' @param method unsigned int
#' @details
#' This is an internal function, it will be excluded in future versions.
#' !!!!! DA IMPLEMENTARE !!!!!
#' @export
sunMC_parBOOT_Cpp <- function(x, B, method) {
    .Call('_INAr_sunMC_parBOOT_Cpp', PACKAGE = 'INAr', x, B, method)
}

#' PIT bootstrap version of the Sun McCabe score test.
#' @param x NumericVector
#' @param B int
#' @param method unsigned int
#' @details
#' This is an internal function, it will be excluded in future versions.
#' !!!!! DA IMPLEMENTARE !!!!!
#' @export
sunMC_pitBOOT_Cpp <- function(x, B, method) {
    .Call('_INAr_sunMC_pitBOOT_Cpp', PACKAGE = 'INAr', x, B, method)
}

sunMCtest_boot <- function(X, arrival, method, B) {
    .Call('_INAr_sunMCtest_boot', PACKAGE = 'INAr', X, arrival, method, B)
}

INARp_cpp <- function(resid, a) {
    .Call('_INAr_INARp_cpp', PACKAGE = 'INAr', resid, a)
}

Xresid <- function(X, alphas, mINN, vINN) {
    .Call('_INAr_Xresid', PACKAGE = 'INAr', X, alphas, mINN, vINN)
}

Xmoments <- function(X, alphas) {
    .Call('_INAr_Xmoments', PACKAGE = 'INAr', X, alphas)
}

