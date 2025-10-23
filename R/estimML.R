
#' Maximum Likelihood for INAR(p) parameter estimation
#'
#' Internal function
#'
#' @param X, observed series
#' @param p, number of lags
#' @param inn, innovation distribution
#' @param control, list of control parameters for optim function
#'
#' @details
#' Descrizione del metodo e Reference alla procedura (Du and Li)
#' @importFrom stats optim
#' @references
#'   \insertAllCited{}
#' @noRd
estimCML <- function(X, p, inn = "poi", control = list()) {
    stopifnot(p==1)
    stopifnot(inn == "poi")

    if (is.null(control$init)){
        YW <- estimYW(X, p, inn = inn)
        alpha0 <- YW$alphas
        lambda0 <- YW$par[1]
        theta0 <- par_transform(alpha0, lambda0, inn)
    }else{
        theta0 <- control$init
    }

    # initial parameter estimation via Yule-Walker

    if(p == 1 & inn == "poi"){
        opt <- optim(theta0, nll_transformed_inar1_poi_ml, data = X, method = "BFGS", control = control) #, hessian = hessian)
        theta_hat <- opt$par
        par_hat <- par_back(theta_hat,inn)

        # # standard errors via Hessian by chatGPT CHECK!!!
        # if (hessian && !is.null(opt$hessian)) {
        #     H <- opt$hessian
        #     if (all(is.finite(H)) && det(H) != 0) {
        #         cov_theta <- solve(H)
        #         a_hat <- par_hat["alpha"]; l_hat <- par_hat["lambda"]
        #         J <- matrix(c(a_hat * (1 - a_hat), 0, 0, l_hat), nrow = 2)
        #         cov_par <- J %*% cov_theta %*% t(J)
        #         se <- sqrt(diag(cov_par))
        #         names(se) <- c("alpha", "lambda")
        #         out$se <- se
        #     } else {
        #         out$se <- c(alpha = NA, lambda = NA)
        #     }
        # }
    }

    OUT <- list("alphas" = par_hat$alphas,
                "par" = par_hat$par
                # "meanX" = mX, "varX" = vX,
                # "meanINN" = mINN, "varINN" = vINN
                )
    return(OUT)
}

# esempio
# set.seed(123)
# x <- genINAR(100000,a = 0.5, par = 2,arrival = "poisson",burnout = 500)$X
# fit <- estimCML(x, p = 1, inn = "poi")
# fit


#' Transformed log-likelihood for MLE
#'
#' @description
#' Transformed negative log-likelihood for INAR(1)-Poisson
#'
#'
#' @keywords internal
#' @noRd
nll_transformed_inar1_poi_ml <- function(data, theta) {
    X <- data
    par_b <- par_back(theta, inn = "poi")
    alpha <- par_b$alphas
    lambda <- par_b$par
    if (alpha <= 0 || alpha >= 1 || lambda <= 0) return(1e10)
    ll <- inar1_poi_loglik_cpp(X, alpha, lambda)
    if (!is.finite(ll)) return(1e10)
    return(-ll)
}

