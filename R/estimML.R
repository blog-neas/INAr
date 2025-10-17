
#' Maximum Likelihood for INAR(p) parameter estimation
#'
#' Internal function
#'
#' @param x, observed series
#' @param p, number of lags
#'
#' @importFrom stats acf
#' @importFrom stats var
#' @importFrom RcppML nnls
#'
#' @details
#' Reference alla procedura (Du and Li)
#' @references
#'   \insertAllCited{}
#' @noRd
estimML <- function(x, p, ...) {
    n <- length(x)
    mX <- mean(x)
    vX <- var(x)

    # ''''''''''''''''''''''''''
    # Some magic happens here
    # ..........................
    OUT <- list(alphas = alphas,
                "meanX" = mX, "varX" = vX,
                "meanINN" = mINN, "varINN" = vINN)
    return(OUT)
}


# versione chatGPT

# =========================================================
# INAR(1)-Poisson: Exact MLE via Rcpp
# =========================================================

library(Rcpp)

# ---------------------------------------------------------
# Parte C++: log-verosimiglianza
# ---------------------------------------------------------
# Rcpp::sourceCpp(code = '
# #include <Rcpp.h>
# using namespace Rcpp;
#
# // [[Rcpp::export]]
# double inar1_poi_loglik_cpp(const IntegerVector& x, double alpha, double lambda) {
#   int n = x.size();
#   if (n < 2) return NA_REAL;
#   double ll = 0.0;
#
#   for (int t = 1; t < n; t++) {
#     int xt = x[t];
#     int xprev = x[t - 1];
#     if (xt < 0 || xprev < 0) return NA_REAL;
#
#     int kmax = std::min(xt, xprev);
#     std::vector<double> log_p(kmax + 1);
#
#     for (int k = 0; k <= kmax; k++) {
#       // combinazione log
#       double lchoose_val = R::lchoose(xprev, k);
#       // log densità Poisson (x_t - k)
#       int diff = xt - k;
#       if (diff < 0) { log_p[k] = -INFINITY; continue; }
#
#       double lp = lchoose_val + k * log(alpha) +
#         (xprev - k) * log(1.0 - alpha) +
#         R::dpois(diff, lambda, true);
#       log_p[k] = lp;
#     }
#
#     // log-sum-exp trick per stabilità numerica
#     double m = -INFINITY;
#     for (int k = 0; k <= kmax; k++) if (log_p[k] > m) m = log_p[k];
#     double sum_exp = 0.0;
#     for (int k = 0; k <= kmax; k++) sum_exp += exp(log_p[k] - m);
#     double log_p_xt = m + log(sum_exp);
#
#     if (!R_finite(log_p_xt)) return -INFINITY;
#     ll += log_p_xt;
#   }
#   return ll;
# }
# ')

# ---------------------------------------------------------
# Wrapper R: ottimizzazione MLE
# ---------------------------------------------------------
estimML <- function(X, p, inn = "poi", control = list()) {
    stopifnot(p==1)
    stopifnot(inn == "poi")

    # if (is.null(init)) {
    #     YW <- estimYW(X, p)
    #     mu <- YW$meanX
    #     alpha0 <- YW$alphas
    #     lambda0 <- max(0.1, mu * (1 - alpha0))
    #     init <- c(alpha = alpha0, lambda = lambda0)
    # # }

    # initial parameter estimation via Yule-Walker
    YW <- estimYW(X, p, inn = inn)
    alpha0 <- YW$alphas
    lambda0 <- YW$par[1]
    theta0 <- par_transform(alpha0, lambda0, inn)

    if(p == 1 & inn == "poi"){
        opt <- optim(theta0, nll_transformed_inar1_poi, method = "BFGS", control = control) #, hessian = hessian)
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
set.seed(123)
x <- genINAR(200,a = 0.5, par = 2,arrival = "poisson",burnout = 500)$X
fit <- estimML(x, p = 1, inn = "poi")

# Parameter trick:
# 1. to_theta: par_transform: transform parameters to unconstrained space
#    alpha in (0,1) -> logit(alpha) = log(alpha/(1-alpha)) in R
#    lambda in (0,inf) -> log(lambda) in R
# 2. from_theta: par_back: transform back to original space
#    logit(alpha) -> alpha = exp(logit(alpha)) / (1 + exp(logit(alpha)))
#    log(lambda) -> lambda = exp(log(lambda))
# how to use:
# 1. transform initial parameters
# 2. transform back the parameters within the optimization function to compute the nll
#    in this way the nll uses the original parameters
#    here the thick is to compute a classical nll function but call it in a nll_transformed
#    which is the one passed to optim and has the only task to transform back the parameters
#    example:
# nll_transformed <- function(theta) {
#     par <- from_theta(theta)  # passa da spazio non vincolato a quello reale
#     alpha <- par["alpha"]
#     lambda <- par["lambda"]
#
#     # valuta la negative log-likelihood
#     saddle_nll_cpp(x, alpha, lambda)
# }
# 3. nll is computed by using real parameters but nll_transformed is linked to
# the transformed parameters, so to obtain the real estimates we need
# to transform back the resulting parameter coming out from optim

# trasformazioni parametri
# to_theta <- function(par) c(qlogis(par["alpha"]), log(par["lambda"]))
# from_theta <- function(theta) c(alpha = plogis(theta[1]), lambda = exp(theta[2]))


#' Transformed log-likelihood
#'
#' @description
#' Transformed negative log-likelihood for INAR(1)-Poisson
#'
#'
#' @keywords internal
#' @noRd
nll_transformed_inar1_poi <- function(theta) {
    par_b <- par_back(theta, inn = "poi")
    alpha <- par_b$alphas
    lambda <- par_b$par
    if (alpha <= 0 || alpha >= 1 || lambda <= 0) return(1e10)
    ll <- inar1_poi_loglik_cpp(x, alpha, lambda)
    if (!is.finite(ll)) return(1e10)
    return(-ll)
}

