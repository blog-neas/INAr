# INARfit
# fit INAR(p) model

#' Fitting INAR(p) Models
#'
#' @param X data vector
#' @param p the order of the INAR(p) process
#' @param inn distribution of the innovation process
#' @param method estimation method
#'
#' @return The fitted model, an object of class `INAR`
#' @details
#' Frontend function that estimates an INAR(p) model given the distribution of the innovation process.
#' @export
INAR <- function(X, p, inn="poi", method = "CLS"){
    # cl <- match.call()
    n <- length(X)
    inn <- trimws(tolower(inn))

    stopifnot(all(X == as.integer(X)))
    stopifnot(p < n)
    if(inn == "negbin" & var(X) <= mean(X)){ stop( "Only overdispersed data allowed for the negbin case" ) }

    if(method == "YW"){
        est <- estimYW(X, p, inn = inn)
    }else if(method == "CLS"){
        est <- estimCLS(X, p, inn = inn)
        alphas <- est$alphas
    }else if(method == "CML"){
        est <- estimCLM(X, p, inn = inn)
        alphas <- est$alphas
        # TO DO
        stop("CML: TO DO")
    }else if(method == "SP"){
        # TO DO
        stop("SP: TO DO")
    }else{ stop('Specify a valid method. Available options: "YW", "CLS", "CML", "SP"') }

    # Innovation Parameters
    innest <- est_pars(est$meanX, est$varX, est$meanINN, est$varINN, inn)

    alphas <- est$alphas
    attr(alphas,"names") <- paste0("a",1:p)

    pars <- innest$pars
    # attr(moments,"pars") # lo fa già est_mom

    momentsINN <- c(est$meanINN, est$varINN)
    attr(momentsINN,"names") <- c("mean", "sigma2")

    coef <- list(
        "alphas"=alphas,
        "pars"=pars
    )

    # fitted lo escludo perchè faccio funzione esterna
    # fitted <- genINAR(n, a = alphas, par = pars, inn = inn, burnout = 500)$X
    resid <- Xresid(X = X, alphas = alphas, mINN = est$meanINN, vINN = est$varINN)
    rms <- sqrt(mean(resid$resid^2,na.rm = TRUE))

    # TO DO
    # - varianza stimatori
    # - loglik e bic

    #AIC <-  n*log(???) + 2*p ????
    #BIC <-  n*log(???) + p*log(n) ????

    # TO DO! --------------------------------
    var.alphas <- matrix(NA,p,p)
    var.pars <- innest$vcov

    var <- list(
        "alphas"=var.alphas,
        "pars"=var.pars
    )

    mask <- list(
        "alphas"=rep(TRUE,p),
        "pars"=rep(TRUE,length(pars))
    )

    resid <- list(
        "res.raw"=resid$resid,
        "res.std"=resid$stdresid
    )

    OUT <- list(
        call = match.call(), call0 = paste0("INAR(",p,") model with ",inn," innovations"),
        data = X, n = n, momentsINN = momentsINN, inn = inn,
        coef = coef, var.coef = var, mask = mask,
        loglik = NA, aic = NA, bic = NA, RMSE = rms,
        residuals = resid, SMCtest = NULL,
        model = paste0(toupper(inn),"-INAR(",p,")")
    )
    class(OUT) <- "INAR" # structure(OUT, class = "INAR")
    return(OUT)
}


#' Estimation of the innovation process' parameters
#'
#' @param mX mean of the INAR process
#' @param varX variance of the INAR process
#' @param mINN mean of the innovation process
#' @param varINN variance of the innovation process
#' @param inn distribution of the innovation process
#'
#' @return A list with parameter estimates
#' @details
#' Inner function that estimates the parameters related with the innovation process, given anINAR(p) model.
#' @noRd
est_pars <- function(mX,varX,mINN,varINN,inn){
    # parameter estimation for CML and YW methods

    if(inn == "poi"){
        lambda <- mINN
        pars <- lambda
        attr(pars,"names") <- "lambda"

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(inn == "negbin"){
        diffvarmu <- abs(varINN - mINN) # trick
        gamma <- (mINN^2)/diffvarmu
        # pi <- mINN/varINN # old
        pi <- diffvarmu/varINN

        pars <- c(gamma, pi)
        attr(pars,"names") <- c("gamma", "pi")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(inn == "genpoi"){
        kappa <- 1 - sqrt(mINN/varINN)
        lambda <- mINN*sqrt(mINN/varINN)

        pars <- c(lambda, kappa)
        attr(pars,"names") <- c("lambda", "kappa")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(inn == "geom"){
        pi <- 1/(1 + mINN)

        pars <- c(pi)
        attr(pars,"names") <- c("pi")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(inn == "yule"){
        rho <- 1/(mINN - 1)

        pars <- c(rho)
        attr(pars,"names") <- c("rho")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(inn == "poislind"){
        # Lívio, T., Khan, N. M., Bourguignon, M., & Bakouch, H. S. (2018).
        # - An INAR (1) model with Poisson–Lindley innovations. Econ. Bull, 38(3), 1505-1513.
        mINN_1 <- 1/mINN
        theta <- (mINN_1 - 1 + sqrt(mINN_1^2 + 6*mINN_1 + 1))/2


        pars <- c(theta)
        attr(pars,"names") <- c("theta")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }

    OUT <- list("pars"=pars,"vcov"=vcov)
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
