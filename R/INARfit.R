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
    par_inn <- getMINN(list(alphas=a_hat,meanX=mean(X),varX=var(X)), inn)

    # resid <- Xresid(X = X, alphas = a_hat, mINN = est$meanINN, vINN = est$varINN)
    # RMSE <- sqrt(mean(resid$resid^2,na.rm = TRUE))
    # mean innovations

    fitted <- INARfitted_cpp(X, par_inn$mINN, a_hat)
    residuals <- X - fitted

    OUT <- list(
        "alphas" = a_hat,
        "par" = par_hat,
        "residuals" = residuals,
        "fitted.values" = fitted
    )
    # class(OUT) <- "INAR" # structure(OUT, class = "INAR")
    return(OUT)
}

#' INAR(p) innovation estimation
#'
#' Internal function
#'
#' @param est list, containing alphas, meanX and varX, output of a estimation procedure (YW, CLS, CML, SP)
#' @param inn character, distribution of the innovation process
#'
#' @details
#' Function that estimates the parameters of the innovation process given
#' the estimated alphas and the moments of the observed series. It is used
#' within the YW procedure to get the innovation parameters, and is also used
#' as initial values for the other estimation procedures (CML, SP).
#'
#' @noRd
getMINN <- function(est, inn){
    # est = output di un metodo di stima (es. YW, CLS, CML, SP)
    stopifnot(inn %in% info_inn$inn)

    alphas <- est$alphas
    mX <- est$meanX
    vX <- est$varX
    R <- est$R

    # innovation moments from X moments and alphas
    # if(is.na(mINN)) mINN <- (1 - sum(alphas))*mX
    # if(is.na(vINN)) vINN <- vX*(1 - sum(alphas^2)) - mX*sum(alphas*(1-alphas))
    mINN <- (1 - sum(alphas))*mX
    if(length(alphas) > 1){
        vINN <- vX*(1 - t(alphas)%*%R%*%alphas) - mX*sum(alphas*(1-alphas))
    }else{
        vINN <- vX*(1 - sum(alphas^2)) - mX*sum(alphas*(1-alphas))
    }

    # OLD
    # if(inn == "poi"){
    #     par <- c("lambda" = mINN)
    # }else if(inn == "negbin"){
    #     diffvarmu <- abs(vINN - mINN)
    #     gamma <- (mINN^2)/diffvarmu
    #     pi <- diffvarmu/vINN
    #
    #     par <- c("gamma" = gamma, "pi" = pi)
    # }else if(inn == "genpoi"){
    #     # CHECK
    #     kappa <- 1 - sqrt(mINN/vINN)
    #     lambda <- mINN*sqrt(mINN/vINN)
    #
    #     par <- c("lambda" = lambda, "kappa" = kappa)
    # }else if(inn == "katz"){
    #     # TO DO
    # }else{
    #     stop("Innovation distribution not implemented yet.")
    # }

    OUT <- list("mINN" = mINN, "vINN" = vINN)
    return(OUT)
}

#' INAR(p) innovation parameters estimation procedures
#'
#' Internal function
#'
#' @param mINN numeric, mean of the innovation process
#' @param vINN numeric, variance of the innovation process
#' @param inn character, distribution of the innovation process
#'
#' @details
#' Function that estimates the parameters of the innovation process given
#' the estimated alphas and the moments of the observed series. It is used
#' within some estimation procedures (YW, CLS).
#'
#' @noRd
getPAR <- function(mINN, vINN, inn){
    if(inn == "poi"){
        par <- c("lambda" = mINN)
    }else if(inn == "negbin"){
        diffvarmu <- abs(vINN - mINN)
        gamma <- (mINN^2)/diffvarmu
        # pi <- diffvarmu/vINN
        # gamma <- mINN^2/abs(vINN - mINN)
        pi <- mINN/vINN

        # p = mu/s2
        # r = mu^2/ (s2 - mu)


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
    return(par)
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
