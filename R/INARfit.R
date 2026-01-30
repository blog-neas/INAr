# INARfit
# fit INAR(p) model

#' Fitting INAR(p) Models
#'
#' @param X data vector
#' @param p the order of the INAR(p) process
#' @param inn distribution of the innovation process, one of "poi" (Poisson), "negbin" (Negative Binomial), "genpoi" (Generalized Poisson), "katz" (Katz)
#' @param method estimation method, one of "YW" (Yule-Walker), "CLS" (Conditional Least Squares), "CML" (Conditional Maximum Likelihood), "SP" (Saddlepoint Approximation)
#'
#' @return The fitted model, an object of class `INAR`
#' @details
#' This function estimates an INAR(p) model given the distribution of the innovation process.
#' @export
INAR <- function(X, p, inn="poi", method = "CLS"){
    # cl <- match.call()
    n <- length(X)
    inn <- trimws(tolower(inn))
    stopifnot(inn %in% info_inn$inn)

    stopifnot(all(X == as.integer(X)))
    stopifnot(p < n)
    if(inn == "negbin" & var(X) <= mean(X)){ stop( "Data are underdispersed. Only overdispersed data allowed for the negbin case" ) }

    if(method == "YW"){
        est <- estimYW(X, p, inn = inn)
    }else if(method == "CLS"){
        est <- estimCLS(X, p, inn = inn)
    }else if(method == "CML"){
        est <- estimCML(X, p, inn = inn)
    }else if(method == "SP"){
        est <- estimSP(X, p, inn = inn)
    }else{
        stop('Specify a valid method. Available options: Yule-Walker "YW",
             Conditional Least Squares "CLS",
             Conditional Maximum Likelihood "CML", Saddlepoint "SP"')
    }

    a_hat <- est$alphas
    par_hat <- est$par

    # resid <- Xresid(X = X, alphas = a_hat, mINN = est$meanINN, vINN = est$varINN)
    # RMSE <- sqrt(mean(resid$resid^2,na.rm = TRUE))

    OUT <- list(
        "alphas" = a_hat,
        "par" = par_hat
    )
    # class(OUT) <- "INAR" # structure(OUT, class = "INAR")
    return(OUT)
}

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
estimPAR <- function(alphas, mX, vX, mINN = NA, vINN = NA, inn = "poi"){
    stopifnot(inn %in% info_inn$inn)

    # innovation moments from X moments and alphas
    if(is.na(mINN)) mINN <- (1 - sum(alphas))*mX
    if(is.na(vINN)) vINN <- vX*(1 - sum(alphas^2)) - mX*sum(alphas*(1-alphas))

    if(inn == "poi"){
        par <- c("lambda" = mINN)
    }else if(inn == "negbin"){
        # CHECK
        diffvarmu <- abs(vINN - mINN) # trick
        gamma <- (mINN^2)/diffvarmu
        # pi <- mINN/vINN # old
        pi <- diffvarmu/vINN

        par <- c("gamma" = gamma, "pi" = pi)
    }else if(inn == "genpoi"){
        # CHECK
        kappa <- 1 - sqrt(mINN/vINN)
        lambda <- mINN*sqrt(mINN/vINN)

        pars <- c("lambda" = lambda, "kappa" = kappa)
    }else if(inn == "katz"){
        # TO DO
    }else{
        stop("Innovation distribution not implemented yet.")
    }

    OUT <- list("par" = par)
    return(OUT)
}


# TENERE SEMPRE COMMENTATO!
# veloce esempio --------------------------------------------------------
# N <- 500
# y <- genINAR(N,0.1,par=1.2,inn="poisson")$X
# INAR(y, p=1)
# y <- genINAR(N,c(0.9,0.01),par=2,inn="poisson")$X
# INARfit(y, p=2)
# y <- genINAR(N,0.1,par=c(1,0.5),inn="negbin")$X
# INARfit(y, p=1, inn="negbin)
# y <- genINAR(N,c(0.9,0.01),par=c(2,0.66),inn="poisson")$X
# INARfit(y, p=2, inn="negbin)
