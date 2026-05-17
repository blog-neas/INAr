#' Yule-Walker for INAR(p) alphas parameter estimation
#'
#' Internal function
#'
#' @param X, observed series
#' @param p, number of lags
#' @param inn, distribution of the innovation process
#' @param ..., additional parameters
#'
#' @importFrom stats acf var
#' @importFrom RcppML nnls
#'
#' @details
#' Reference alla procedura (Du and Li)
#' @references
#'   \insertAllCited{}
#' @noRd
estimYW <- function(X, p, inn = "poi", ...) {
    stopifnot(inn %in% info_inn$inn)

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
        R <- 1
    }

    names(alphas) <- paste0("a",1:p)

    INN_par <- getMINN(list(alphas=alphas,meanX=mean(X),varX=var(X), R=R),inn)
    par <- getPAR(INN_par$mINN, INN_par$vINN, inn)


    OUT <- list("alphas" = alphas,
                "par"=par,
                "meanX" =  mean(X), "varX" = var(X)
    )
    return(OUT)
}


# library(INAr)
# ii <- "negbin"
# pp <- c(5, 0.8) # per negbin
# aa <- genINAR(100000, a = 0.4, par = pp, inn = ii)$X
# INAr:::estimYW(aa, p = 1, inn = ii)$par;pp
# INAr:::estimCLS(aa, p = 1, inn = ii)$par;pp
# bb <- genINAR(100000, a = c(0.4,0.2), par = pp, inn = ii)$X
# INAr:::estimYW(bb, p = 2, inn = ii)$par;pp
# INAr:::estimCLS(bb, p = 2, inn = ii)$par;pp
# cc <- genINAR(100000, a = c(0.4,0.2,0.15), par = pp, inn = ii)$X
# INAr:::estimYW(cc, p = 3, inn = ii)$par;pp
# INAr:::estimCLS(cc, p = 3, inn = ii)$par;pp


# library(INAr)
# xx <- genINAR(100000, a = 0.5, par = 2, inn = "poi")$X
# INAr:::estimYW(xx, p = 1, inn = "poi")
# zz <- genINAR(100000, a = 0.4, par = c(5,0.8), inn = "negbin")$X
# INAr:::estimYW(zz, p = 1, inn = "negbin")
# yy <- genINAR(100000, a = 0.3, par = c(2,0.5), inn = "genpoi")$X
# INAr:::estimYW(yy, p = 1, inn = "genpoi")
# #
# xx <- genINAR(100000, a = c(0.5,0.2), par = 2, inn = "poi")$X
# INAr:::estimYW(xx, p = 2, inn = "poi")
# zz <- genINAR(100000, a = c(0.5,0.2), par = c(5,0.8), inn = "negbin")$X
# INAr:::estimYW(zz, p = 2, inn = "negbin")
# yy <- genINAR(100000, a = c(0.5,0.1), par = c(2,0.5), inn = "genpoi")$X
# INAr:::estimYW(yy, p = 2, inn = "genpoi")

# library(INAr)
# xx <- genINAR(100000, a = 0.5, par = 2, inn = "poi")$X
# INAr:::estimCLS(xx, p = 1, inn = "poi")
# zz <- genINAR(100000, a = 0.4, par = c(5,0.8), inn = "negbin")$X
# INAr:::estimCLS(zz, p = 1, inn = "negbin")
# yy <- genINAR(100000, a = 0.3, par = c(2,0.5), inn = "genpoi")$X
# INAr:::estimCLS(yy, p = 1, inn = "genpoi")
# #
# xx <- genINAR(100000, a = c(0.5,0.2), par = 2, inn = "poi")$X
# INAr:::estimCLS(xx, p = 2, inn = "poi")
# zz <- genINAR(100000, a = c(0.5,0.2), par = c(5,0.8), inn = "negbin")$X
# INAr:::estimCLS(zz, p = 2, inn = "negbin")
# yy <- genINAR(100000, a = c(0.5,0.1), par = c(2,0.5), inn = "genpoi")$X
# INAr:::estimCLS(yy, p = 2, inn = "genpoi")

