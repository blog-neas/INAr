# INARfit
# fit INAR(p) model

INARfit <- function(X,order,arrival="poisson"){

    n <- length(X)
    err <- NULL

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
            names(err) <- paste0("fcnnls_solve_",1:order)
        }
    }else{
        a <- r
    }

    # sempre secondo Du and Li
    resid <- rep(NA,n-order)
    for(t in (order+1):n ){
        resid[t - order] <- X[t] - X[(t-1):(t-order)]%*%a
    }
    mX <- mean(X)
    varX <- var(X)
    mINN <- (1 - sum(a))*mean(X)
    varINN <- (1/(n-order))*sum((resid - mean(resid))^2)

    arr_mom <- est_mom(mX,varX,mINN,varINN,arrival)

    names(a) <- paste0("a",1:order)
    out <- list("alphas"=a,"par"=arr_mom,"moments"=c("mean_X"=mX,"var_X"=varX,"mean_EPS"=mINN,"var_EPS"=varINN),"resid"=resid,"err"=err)
    return(out)
}

est_mom <- function(mX,varX,mINN,varINN,arrival){
    # if(arrival=="poisson"){}
    lambda <- mINN
    return(c("lambda"=lambda))
}

# veloce esempio --------------------------------------------------------
# N <- 500
# y <- genINAR(N,0.1,0.3,arrival="poisson")$X
# INARfit(y,order=1)
# y <- genINAR(N,c(0.9,0.01),2,arrival="poisson")$X
# INARfit(y,order=2)
