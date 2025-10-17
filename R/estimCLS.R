#' Conditional Nonnegative Least Squares for INAR(p) parameter estimation
#'
#' Internal function
#'
#' @param x vector, observed series
#' @param p integer, number of lags
#'
#' @importFrom stats var
#' @importFrom RcppML nnls
#'
#' @details
#' Reference alla procedura
#' @references
#'   \insertAllCited{}
#' @noRd
estimCLS <- function(x, p) {
    n <- length(x)
    mX <- mean(x)
    vX <- var(x)

    Yreg <- x[(p+1):n]
    # Xreg <- lagmat(x,p)

    Xreg <- lagmat(x,p)

    ## nnls::nnls
    # mod <- nnls(Xreg, Yreg)
    # alphas <- mod$x[-1]
    # attr(alphas, "names") <- paste0("a",1:p)
    #
    # mINN <- mod$x[1]
    # vINN <- mod$deviance/(n-(p+1))

    ## RcppML::nnls
    mod <- nnls(crossprod(Xreg), crossprod(Xreg, Yreg))
    alphas <- mod[-1]
    attr(alphas, "names") <- paste0("a",1:p)

    mINN <- mod[1]
    vINN <- sum((Yreg - Xreg%*%mod)^2)/(n-(p+1))

    par_hat <- estimPAR(alphas, mX, vX, mINN, vINN, inn = "poi")

    # est <- est_mom(arr_mom$meanX, arr_mom$varX, arr_mom$meanINN, arr_mom$varINN, arrival)

    OUT <- list(alphas = alphas,
                par = par_hat$par,
                "meanX" = mX, "varX" = vX,
                "meanINN" = mINN, "varINN" = vINN)
    return(OUT)
}

