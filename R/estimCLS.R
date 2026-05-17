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
estimCLS <- function(X, p, inn){
    n <- length(X)

    Yreg <- X[(p+1):n]
    Xreg <- lagmat(X,p)

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

    # OLD, MA ERRATO, VARIANCE INNOVAZIONI != VARIANCE RESIDUI
    vINN <- sum((Yreg - Xreg%*%mod)^2)/(n-(p+1))


    # residui
    res <- as.numeric(Yreg - Xreg %*% mod)

    # var dei residui != var. innovazioni
    res_var <- mean(res^2)

    # medie dei thinning laggati
    lag_means <- colMeans(Xreg[, -1, drop = FALSE])

    thinning_var <- sum(alphas * (1 - alphas) * lag_means)

    vINN <- res_var - thinning_var

    if(!is.finite(vINN) || vINN <= 0){
        stop("Estimated innovation variance must be finite and strictly positive.")
    }
    if(inn == "negbin" && (!is.finite(mINN) || vINN <= mINN)){
        stop("Negative-binomial innovations require estimated variance to be finite and strictly greater than the estimated mean.")
    }

    par <- getPAR(mINN, vINN, inn)

    # if(inn == "poi"){
    #     par <- c("lambda" = mINN)
    # }else if(inn == "negbin"){
    #     # CHECK!
    #     diffvarmu <- abs(vINN - mINN)
    #     gamma <- (mINN^2)/diffvarmu
    #     pi <- diffvarmu/vINN
    #
    #     par <- c("gamma" = gamma, "pi" = pi)
    # }else if(inn == "genpoi"){
    #     kappa <- 1 - sqrt(mINN/vINN)
    #     lambda <- mINN*sqrt(mINN/vINN)
    #
    #     par <- c("lambda" = lambda, "kappa" = kappa)
    # }else if(inn == "katz"){
    #     # TO DO
    # }else{
    #     stop("Innovation distribution not implemented yet.")
    # }

    OUT <- list("alphas" = alphas,
                "par" = par,
                "meanX" = mean(X), "varX" = var(X)
                # "meanINN" = mINN, "varINN" = vINN
                )
    return(OUT)
}



# esempio, tipo unit root test
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

