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


