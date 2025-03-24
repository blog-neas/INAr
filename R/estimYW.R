
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
estimYW <- function(x, p, ...) {
    n <- length(x)
    err <- NULL
    stopifnot(p < n)

    r <- acf(x, plot = FALSE)$acf[2:(p+1)]
    if(p > 1){
        # OLD, sostituito da YW_cpp
        # R <- diag(p)
        # for(i in 1:(p-1)){
        #     for(j in (i+1):p){
        #         R[i,j] <- r[abs(i-j)]
        #     }
        # }
        # R[lower.tri(R)] <- t(R)[lower.tri(R)]

        R <- YW_cpp(r)

        # ALPHAS = R^-1r
        # alphas <- solve(R, r)

        # versione fast nonnegative factorization
        # RcppML::nnls
        alphas <- as.vector(nnls(as.matrix(R),as.matrix(r)))

        # nnls::nnls
        # alphas <- nnls(R, r)$x

        attr(alphas, "names") <- paste0("a",1:p)

        # if(all(alphas > 0)){
        #     err <- NULL
        # }else{
        #     warning("Y-W estimation not reliable. Please rely on more consistent estimators.")
        #     # NMF::fcnnls
        #     # Wtmp[,j] <- RcppML::nnls(crossprod(Xdum), crossprod(Xdum, Adum[,j]))
        #     # err <- as.vector(fcnnls(R,r)$x) # ,pseudo=TRUE per Monroe-Penrose version
        #     names(err) <- paste0("fcnnls_solve_a",1:p)
        # }
    }else{
        alphas <- r
    }

    mX <- mean(x)
    vX <- var(x)
    mINN <- (1 - sum(alphas))*mX
    vINN <- vX*(1 - sum(alphas^2)) - mX*sum(alphas*(1-alphas))
    #  acf(x, plot = FALSE)$acf[1] - sum(alphas*r)

    OUT <- list("alphas" = alphas,
                "meanX" = mX, "varX" = vX,
                "meanINN" = mINN, "varINN" = vINN)
    return(OUT)
}

