#' Perform INAR tests
#'
#' @param X vector, thinning parameters. The lenght of this vector defines the number of lags `p` of the INAR(p) process.
#' @param type character, the type of test to be performed. The options are "smc", "hmc" and "rho".
#' @param B integer, the number of bootstrap samples to be generated. The default value is 0, corresponding to no bootstrap samples.
#' @return A number.
#'
#' @details
#' The function performs tests for the INAR(1) model.
#' The tests are based on the Sun-McCabe (SMC), Harris-McCabe (HMC) and Rho statistics.
#' The tests can be performed using the original data or using bootstrap samples.
#' The bootstrap samples are generated using the parametric, semiparametric and nonparametric methods.
#' It is possible to perform the tests using the original data (B=0) or using bootstrap samples (B > 0).
#' @examples
#' # ....... examples .....
#' # INARtest(X = rpois(100,2), type = "smc", B = 0)
#'
#' @export
INARtest <- function(X, type, B = 0){
    type <- tolower(type)
    stopifnot(type %in% c("smc","hmc","rho"))
    OUT <- NA # completa lo script

    if(B > 1){
        # perform bootstrap tests
        if(type == "smc"){
            # call smc.test
        }else if(type == "hmc"){
            # call hmc.test
        }else if(type == "rho"){
            # call rho.test
        }

    }else{
        # perform simple tests
        if(type == "smc"){
            # call smc.test
        }else if(type == "hmc"){
            # call hmc.test
        }else if(type == "rho"){
            # call rho.test
        }

    }

    # TO DO:
    # configure OUT to have a nice output in the style of a summary of a test

    return(OUT)
}

#' Perform Sun-McCabe INAR tests
#'
#' @param X vector, thinning parameters. The lenght of this vector defines the number of lags `p` of the INAR(p) process.
#' @param method character, the distribution to be used in the test. The options are "poi", "negbin" and "genpoi".
#' @param type character, the type of test to be performed, the alternatives are "parametric or "semiparametric".
#' @param B integer, the number of bootstrap samples to be generated. The default value is 0, corresponding to no bootstrap samples.
#' @return A number.
#'
#' @details
#' The function performs the Sun-McCabe (SMC) tests for the INAR(1) model.
#' It is possible to run an exact (B = 0) or a bootstrap (B > 0) test.
#' SMC tests are based on Poisson, negative binomial and generalized Poisson
#' distributions, in all the considered cases, both parametric and semiparametric
#' methods are available.
#' It is possible to perform the tests using the original data (B=0) or using
#' bootstrap samples (B > 0).
#' @examples
#' # ....... examples .....
#' # SMCtest(X = rpois(100,2), type = "semiparametric", B = 0)
#'
#' @export
SMCtest <- function(X, method, type = NA, B = 0){
    type <- tolower(type)
    stopifnot(method %in% 1:3)

    B_stat <- NA
    B_pval <- NA

    # perform exact test
    smc_est <- SMC_Cpp(X, method)
    stat <- smc_est$stat
    pval <- smc_est$pval

    if(B > 0){
        stopifnot(type %in% c("parametric","semiparametric"))
        # perform bootstrap test
        if(type == "parametric"){
            # call smc.test
            smc_boot <- SMC_parBOOT_Cpp(X, B, method)
        }else if(type == "semiparametric"){
            # call smc.test
            smc_boot <- SMC_semiparBOOT_Cpp(X, B, method)
        }else{
            stop("Specify a correct test type")
        }
        B_stat <- mean(smc_boot)
        B_pval <- mean(abs(smc_boot) > abs(stat), na.rm = TRUE) # imbroglio, uso: na.rm = TRUE
    }

    OUT <- list(
        stat = stat,
        pval = pval,
        B_stat = B_stat,
        B_pval = B_pval
        )

    # TO DO:
    # configure OUT to have a nice output in the style of a summary of a test

    return(OUT)
}


