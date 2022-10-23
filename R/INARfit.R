# INARfit
# fit INAR(p) model

#' Estimation of the innovation process' parameters
#'
#' @param X data vector
#' @param order p, the order of the INAR(p) process
#' @param arrival distribution of the innovation process
#'
#' @return The fitted model
#' @details
#' Inner function that estimates the parameters related with the innovation process, given anINAR(p) model.
#' @export
INARfit <- function(X,order,arrival="poisson"){

    n <- length(X)
    err <- NULL
    stopifnot(order < n)

    # YULE-WALKER ESTIMATION ----------------------------
    # - secondo Du and Li -------------------------------

    # FIT ALPHAS
    # da teoria
    # R %*% ALPHAS = r
    r <- acf(X, plot = FALSE)$acf[2:(order+1)]
    if(order > 1){
        R <- diag(order)
        for(i in 1:(order-1)){
            for(j in (i+1):order){
                R[i,j] <- r[abs(i-j)]
            }
        }
        R[lower.tri(R)] <- t(R)[lower.tri(R)]

        # ALPHAS = R^-1r
        a <- solve(R, r)
        if(all(a>0)){
            err <- NULL
        }else{
            warning("Y-W estimation not reliable. Please rely on more consistent estimators.")
            err <- as.vector(NMF::fcnnls(R,r)$x) # ,pseudo=TRUE per Monroe-Penrose version
            names(err) <- paste0("fcnnls_solve_a",1:order)
        }
    }else{
        a <- r
    }

    # sempre secondo Du and Li
    # resid <- rep(NA,n-order)
    # for(t in (order+1):n){
    #     resid[t - order] <- X[t] - X[(t-1):(t-order)]%*%a
    # }

    # mX <- mean(X)
    # varX <- var(X)
    # mINN <- (1 - sum(a))*mean(X)
    # varINN <- (1/(n-order))*sum((resid - mean(resid))^2)

    arr_mom <- Xmoments(X,a)
    est <- est_mom(arr_mom$meanX, arr_mom$varX, arr_mom$meanINN, arr_mom$varINN, arrival)

    names(a) <- paste0("a",1:order)
    out <- list("alphas"=a,"par"=est,
                "moments"=c("mean_X"=arr_mom$meanX,"var_X"=arr_mom$varX,"mean_EPS"=arr_mom$meanINN,"var_EPS"=arr_mom$varINN),"resid"=arr_mom$resid,"stdresid"=arr_mom$stdresid,"err"=err)
    return(out)
}

#' Estimation of the innovation process' parameters
#'
#' @param mX mean of the INAR process
#' @param vX variance of the INAR process
#' @param mINN mean of the innovation process
#' @param vINN variance of the innovation process
#' @param arrival distribution of the innovation process
#'
#' @return A list with parameter estimates
#' @details
#' Inner function that estimates the parameters related with the innovation process, given anINAR(p) model.
#' @export
est_mom <- function(mX,varX,mINN,varINN,arrival){
    if(arrival == "poisson"){
        lambda <- mINN
        OUT <- list("lambda"=lambda)
    }
    if(arrival == "negbin"){

        # trick
        diffvarmu <- abs(varINN - mINN)
        gamma <- (mINN^2)/diffvarmu
        # pi <- mINN/varINN # old
        pi <- diffvarmu/varINN

        OUT <- list("gamma"=gamma,"pi"=pi)
    }
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
