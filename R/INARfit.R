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
    stopifnot(all(X >= 0))
    stopifnot(p < n)
    if(inn == "negbin" & var(X) <= mean(X)){ stop( "Data are underdispersed. Only overdispersed data allowed for the Negative Binomial case.", "Consider consider using the katz distribution." ) }

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
    par_inn <- getMINN(list(alphas=a_hat,meanX=mean(X),varX=var(X), R = est$R), inn)

    # resid <- Xresid(X = X, alphas = a_hat, mINN = est$meanINN, vINN = est$varINN)
    # RMSE <- sqrt(mean(resid$resid^2,na.rm = TRUE))
    # mean innovations

    fitted <- INARfitted_cpp(X, par_inn$mINN, a_hat)
    residuals <- X - fitted
    # CHECK
    stdresiduals <- residuals/sqrt(as.numeric(par_inn$vINN))


    OUT <- list(
        "call" = match.call(),
        "alphas" = a_hat,
        "par" = par_hat,
        "residuals" = residuals,
        "stdresiduals" = stdresiduals,
        "fitted.values" = fitted,
        "data" = X,
        "inn" = inn,
        "mINN" = par_inn$mINN,
        "vINN" = par_inn$vINN
    )
    class(OUT) <- "INAR" # structure(OUT, class = "INAR")
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
#' @param eps numeric, small positive value to avoid numerical issues
#'
#' @details
#' Function that estimates the parameters of the innovation process given
#' the estimated alphas and the moments of the observed series. It is used
#' within some estimation procedures (YW, CLS).
#'
#' Parameter estimates are obtained from scientific papers and references, according to the respective parametrization of the innovation distribution used in this package. The main references are:
#' - Poisson, see ...
#' - Negative Binomial, see  ...
#' - Generalized Poisson, see  ...
#' - Katz, see  ...
#' - Zero-Inflated Poisson \insertCite{piancastelli2019inferential}{INAr}.
#' - Zero-Inflated Negative Binomial, see  ...
#' - Double Poisson, see ...
#' - Geometric, see ...
#' - Binomial, see ...
#' @references
#'   \insertAllCited{}
#' @return A vector of parameter estimates.
#' @noRd
getPAR <- function(mINN, vINN, inn, eps = 1e-8) {
    if(inn == "poi"){
        par <- c("lambda" = mINN)
    }else if(inn == "negbin") {
        diffvarmu <- abs(vINN - mINN)
        gamma <- (mINN^2)/diffvarmu
        # pi <- diffvarmu/vINN
        # gamma <- mINN^2/abs(vINN - mINN)
        pi <- mINN/vINN

        # p = mu/s2
        # r = mu^2/ (s2 - mu)


        par <- c("gamma" = gamma, "pi" = pi)
    }else if(inn == "genpoi") {
        kappa <- 1 - sqrt(mINN/vINN)
        lambda <- mINN*sqrt(mINN/vINN)

        par <- c("lambda" = lambda, "kappa" = kappa)
    }else if(inn == "katz") {
        #
        if (mINN == 0) {
            return(c(a = 0, b = 0))
        }
        if (vINN == 0) {
            stop("A non-degenerate Katz distribution cannot have positive mean and zero variance.")
        }

        a <- mINN^2 / vINN
        b <- 1 - mINN / vINN

        if (b < 0) {
            m_raw <- -a / b
            m <- round(m_raw)

            if (m < 1 || abs(m_raw - m) > tol) {
                # Projection to the closest binomial-type Katz distribution.
                #
                # For Bin(m, p):
                # mean = m p
                # var  = m p (1 - p)
                #
                # Katz parameters:
                # a = m p / (1 - p)
                # b = -p / (1 - p)
                #
                # From moments:
                # var / mean = 1 - p
                # p = 1 - var / mean
                #
                # m = mean / p

                p <- 1 - vINN / mINN
                p <- min(max(p, tol), 1 - tol)

                m <- max(1L, round(mINN / p))

                p <- min(max(mINN / m, tol), 1 - tol)

                a <- m * p / (1 - p)
                b <- -p / (1 - p)
            }
        }

        par <- c("a" = a, "b" = b)
    }else if(inn == "zip"){
        cvINN <- vINN/mINN
        lambda <- cvINN + mINN - 1
        # sigma = mixing parameter
        sigma <- 1 - mINN/lambda

        par <- c("lambda" = lambda, "sigma" = sigma)
    }else if(inn == "zinb"){
        # TO DO

    }else if(inn == "dpoi"){
        mu <- max(mINN, eps)
        sigma <- max(vINN / mINN, eps)

        par <- c("mu" = mu, "sigma" = sigma)
    }else if (inn == "geom") {
        # Geometrica con supporto {0,1,2,...}: mean = (1-pi)/pi
        # Allora: pi = 1 / (1 + mean)
        pi <- 1 / (1 + mINN)

        par <- c("pi" = pi)
    }else if (inn == "bin") {
        # binomiale con supporto {0,1,...,n}: mean = n*pi
        # Binomial(size, prob): mean = size*prob, variance = size*prob*(1-prob)
        # From moments: prob = 1 - variance/mean, size = mean/prob.

        if (vINN < 0 || vINN > mINN) {
            stop("Binomial innovations require overdispersed data: estimated variance between 0 and the estimated mean.")
        }

        # For rbinom() size must be an integer.
        # Moment estimates will generally not give an exact integer.
        pi <- 1 - vINN / mINN
        enne <- mINN / pi

        # CORREZIONI:
        # enne_adj <- max(1L, as.integer(round(mINN / pi)))
        # # Refit prob after rounding size, keeping mean approximately equal.
        # pi_adj <- min(max(mINN / enne_adj, 0), 1)

        par <- c("n" = enne, "p" = pi)
    }else{
        stop("Innovation distribution not implemented yet.")
    }
    return(par)
}
