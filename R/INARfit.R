# INARfit
# fit INAR(p) model

#' Fitting INAR(p) Models
#'
#' @param X data vector
#' @param order p, the order of the INAR(p) process
#' @param arrival distribution of the innovation process
#' @param method estimation method
#'
#' @return The fitted model, an object of class `INAR`
#' @details
#' Frontend function that estimates an INAR(p) model given the distribution of the innovation process.
#' @export
INAR <- function(X, order, arrival="poisson", method = "CLS"){
    # cl <- match.call()
    n <- length(X)
    arrival <- trimws(tolower(arrival))

    stopifnot(all(X == as.integer(X)))
    stopifnot(order < n)
    if(arrival == "negbin" & var(X) <= mean(X)){ stop( "Only overdispersed data allowed for the negbin case" ) }

    if(method == "YW"){
        mom <- estimYW(X, order)
    }else if(method == "CLS"){
        mom <- estimCLS(X, order)
    }else if(method == "CML"){
        # TO DO
        stop("CML: TO DO")
    }else if(method == "SP"){
        # TO DO
        stop("SP: TO DO")
    }else{ stop('Specify a valid method. Available options: "YW", "CLS", "CML", "SP"') }

    fit_res <- Xresid(X = X, alphas = mom$alphas, mINN = mom$meanINN, vINN = mom$varINN)

    # TO DO
    # - varianza stimatori
    # - loglik e bic

    # Innovation Parameters
    est <- est_pars(mom$meanX, mom$varX, mom$meanINN, mom$varINN, arrival)

    alphas <- mom$alphas
    attr(alphas,"names") <- paste0("a",1:order)

    pars <- est$pars
    # attr(moments,"pars") # lo fa già est_mom

    momentsINN <- c(mom$meanINN, mom$varINN)
    attr(momentsINN,"names") <- c("mean", "sigma2")

    coef <- list(
        "alphas"=alphas,
        "pars"=pars
    )

    # TO DO!
    var.alphas <- matrix(NA,order,order)
    var.pars <- est$vcov

    var <- list(
        "alphas"=var.alphas,
        "pars"=var.pars
    )

    mask <- list(
        "alphas"=rep(TRUE,order),
        "pars"=rep(TRUE,length(pars))
    )

    resid <- list(
        "res.raw"=fit_res$resid,
        "res.std"=fit_res$stdresid
    )

    OUT <- list(
        call = match.call(), call0 = paste0("INAR(",order,") model with ",arrival," innovations"), # call = match.call(),
        data = X, momentsINN = momentsINN, arrival = arrival,
        coef = coef, var.coef = var, mask = mask,
        loglik = NA, aic = NA, bic = NA,
        residuals = resid, SMCtest = NULL,
        model = paste0(toupper(arrival),"-INAR(",order,")")
    )
    class(OUT) <- "INAR" # structure(OUT, class = "INAR")
    OUT
}


#' Estimation of the innovation process' parameters
#'
#' @param mX mean of the INAR process
#' @param varX variance of the INAR process
#' @param mINN mean of the innovation process
#' @param varINN variance of the innovation process
#' @param arrival distribution of the innovation process
#'
#' @return A list with parameter estimates
#' @details
#' Inner function that estimates the parameters related with the innovation process, given anINAR(p) model.
#' @noRd
est_pars <- function(mX,varX,mINN,varINN,arrival){
    # parameter estimation for CML and YW methods

    if(arrival == "poisson"){
        lambda <- mINN
        pars <- lambda
        attr(pars,"names") <- "lambda"

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(arrival == "negbin"){
        diffvarmu <- abs(varINN - mINN) # trick
        gamma <- (mINN^2)/diffvarmu
        # pi <- mINN/varINN # old
        pi <- diffvarmu/varINN

        pars <- c(gamma, pi)
        attr(pars,"names") <- c("gamma", "pi")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(arrival == "genpoi"){
        kappa <- 1 - sqrt(mINN/varINN)
        lambda <- mINN*sqrt(mINN/varINN)

        pars <- c(lambda, kappa)
        attr(pars,"names") <- c("lambda", "kappa")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(arrival == "geom"){
        pi <- 1/(1 + mINN)

        pars <- c(pi)
        attr(pars,"names") <- c("pi")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }else if(arrival == "yule"){
        rho <- 1/(mINN - 1)

        pars <- c(rho)
        attr(pars,"names") <- c("rho")

        # TO DO
        vcov <- matrix(NA,length(pars),length(pars))
    }

    OUT <- list("pars"=pars,"vcov"=vcov)
    return(OUT)
}


# TENERE SEMPRE COMMENTATO!
# veloce esempio --------------------------------------------------------
# N <- 500
# y <- genINAR(N,0.1,par=1.2,arrival="poisson")$X
# INARfit(y, order=1)
# y <- genINAR(N,c(0.9,0.01),par=2,arrival="poisson")$X
# INARfit(y, order=2)
# y <- genINAR(N,0.1,par=c(1,0.5),arrival="negbin")$X
# INARfit(y, order=1, arrival="negbin)
# y <- genINAR(N,c(0.9,0.01),par=c(2,0.66),arrival="poisson")$X
# INARfit(y, order=2, arrival="negbin)
