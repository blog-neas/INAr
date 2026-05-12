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
    vINN <- sum((Yreg - Xreg%*%mod)^2)/(n-(p+1))

    if(inn == "poi"){
        par <- c("lambda" = mINN)
    }else if(inn == "negbin"){
        # CHECK!
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
                "par" = par_hat$par,
                "meanX" = mean(X), "varX" = var(X)
                # "meanINN" = mINN, "varINN" = vINN
                )
    return(OUT)
}
