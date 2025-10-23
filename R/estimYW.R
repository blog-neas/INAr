#' Yule-Walker for INAR(p) alphas parameter estimation
#'
#' Internal function
#'
#' @param x, observed series
#' @param p, number of lags
#'
#' @importFrom stats acf
#' @importFrom stats var
#' @importFrom RcppML nnls
#'
#' @details
#' Reference alla procedura (Du and Li)
#' @references
#'   \insertAllCited{}
#' @noRd
estimYW <- function(X, p, inn = "poi", ...) {
    stopifnot(inn %in% c("poi"))

    n <- length(X)
    err <- NULL
    stopifnot(n > 1, p < n)

    r <- acf(X, plot = FALSE)$acf[2:(p+1)]
    if(p > 1){
        R <- YW_cpp(r)

        # versione fast nonnegative factorization
        # RcppML::nnls
        alphas <- as.vector(nnls(as.matrix(R),as.matrix(r)))

        attr(alphas, "names") <- paste0("a",1:p)
    }else{
        alphas <- r
    }

    names(alphas) <- paste0("a",1:p)

    par <- estimPAR(alphas, mean(X), var(X), inn = inn)

    OUT <- list("alphas" = alphas,
                "par"=par$par
                # "meanX" = mX, "varX" = vX,
                # "meanINN" = mINN, "varINN" = vINN
                )
    return(OUT)
}

