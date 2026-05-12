#' Yule-Walker for INAR(p) alphas parameter estimation
#'
#' Internal function
#'
#' @param X, observed series
#' @param p, number of lags
#' @param inn, distribution of the innovation process
#' @param ..., additional parameters
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

    mINN <- (1 - sum(alphas))*mean(X)
    vINN <- var(X)*(1 - sum(alphas^2)) - mean(X)*sum(alphas*(1-alphas))

    if(inn == "poi"){
        par <- c("lambda" = mINN)
    }else if(inn == "negbin"){
        diffvarmu <- abs(vINN - mINN)
        gamma <- (mINN^2)/diffvarmu
        pi <- diffvarmu/vINN

        par <- c("gamma" = gamma, "pi" = pi)
    }else if(inn == "genpoi"){
        kappa <- 1 - sqrt(mINN/vINN)
        lambda <- mINN*sqrt(mINN/vINN)

        par <- c("lambda" = lambda, "kappa" = kappa)
    }else if(inn == "katz"){
        # TO DO
    }else{
        stop("Innovation distribution not implemented yet.")
    }

    OUT <- list("alphas" = alphas,
                "par"=par,
                "meanX" =  mean(X), "varX" = var(X)
    )
    return(OUT)
}


# generiamo un esempio, tipo unit root test
library(INAR)
xx <- genINAR(1000, a = 0.5, par = 2, inn = "poisson")$X
estimYW(xx, p = 1)


