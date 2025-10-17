#' INAR(p) parameter estimation procedures
#'
#' Internal function
#'
#' @param alphas vector, estimated alphas
#' @param mX numeric, estimated mean of the observed series
#' @param vX numeric, estimated variance of the observed series
#' @param mINN numeric, estimated mean of the innovation process
#' @param vINN numeric, estimated variance of the innovation process
#' @param inn character, distribution of the innovation process
#'
#' @details
#' Function that estimates the parameters of the innovation process given
#' the estimated alphas and the moments of the observed series. It is used
#' within some estimation procedures (YW, CLS).
#'
#' @noRd
estimPAR <- function(alphas, mX, vX, mINN, vINN, inn = "poi"){
    stopifnot(inn %in% c("poi"))

    if(inn == "poi"){
        par = c("lambda" = mINN)
    }
    # continua con altri innovation dgp

    OUT <- list("par"=par)
    return(OUT)
}

