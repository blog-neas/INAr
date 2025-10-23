#' Saddlepoint method for INAR(p) parameter estimation
#'
#' Internal function
#'
#' @param x, observed series
#' @param p, number of lags
#'
#' @importFrom stats var
#' @importFrom RcppML nnls
#'
#' @details
#' Reference alla procedura
#' @importFrom stats optim
#' @references
#'   \insertAllCited{}
#' @noRd
estimSP <- function(X, p, inn = "poi", control = list()){
    stopifnot(p==1)
    stopifnot(inn == "poi")

    if (is.null(control$init)) {
        YW <- estimYW(X, p, inn)
        theta0 <- par_transform(a = YW$alphas,par = YW$par,inn = inn)
    }else{
        theta0 <- par_transform(a = control$init[1:p],par = control$init[-c(1:p)], inn = inn)
    }

    if(p == 1 & inn == "poi"){
        opt <- optim(theta0, nll_transformed_inar1_poi_sp, data = X, method = "BFGS", control = control) #, hessian = hessian)
        theta_hat <- opt$par
        par_hat <- par_back(theta_hat,inn)
    }

    OUT <- list("alphas" = par_hat$alphas,
                "par" = par_hat$par
                # "meanX" = mX, "varX" = vX,
                # "meanINN" = mINN, "varINN" = vINN
    )
    return(OUT)
}


#' Transformed log-likelihood for SPE
#'
#' @description
#' Transformed negative log-likelihood for INAR(1)-Poisson
#'
#'
#' @keywords internal
#' @noRd
nll_transformed_inar1_poi_sp <- function(data, theta) {
    X <- data
    par_b <- par_back(theta, inn = "poi")
    alpha <- par_b$alphas
    lambda <- par_b$par
    if (alpha <= 0 || alpha >= 1 || lambda <= 0 || lambda > 1000) return(1e10)
    ll <- saddle_nll_cpp(X, alpha, lambda)
    if (!is.finite(ll)) return(1e10)
    return(-ll)
}

