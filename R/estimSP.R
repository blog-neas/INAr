
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
#' @references
#'   \insertAllCited{}
#' @noRd
estimSP <- function(X, p) {
    n <- length(X)
    mX <- mean(X)
    vX <- var(X)

    # ''''''''''''''''''''''''''
    # Some magic happens here
    # ..........................

    OUT <- list(alphas = alphas,
                "meanX" = mX, "varX" = vX,
                "meanINN" = mINN, "varINN" = vINN)
    return(OUT)
}

# VERSIONE chatGPT

# ---------------------------------------------------------
# R side
# ---------------------------------------------------------
library(Rcpp)

# Compila automaticamente la parte C++
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double saddle_nll_cpp(const IntegerVector& x, double alpha, double lambda) {
  int n = x.size();
  if (n < 2) return NA_REAL;
  double ll = 0.0;

  for (int t = 1; t < n; t++) {
    int xt = x[t];
    int xprev = x[t - 1];

    if (xt == 0) {
      // exact formula: f_eps(0) * (1 - alpha)^{x_{t-1}}
      double logp = -lambda + xprev * log(1.0 - alpha);
      ll += logp;
      continue;
    }

    double a = lambda * alpha;
    double b = (xprev - xt) * alpha + lambda * (1.0 - alpha);
    double c = - (double)xt * (1.0 - alpha);
    double D = b * b - 4.0 * a * c;
    if (D < 0.0) return 1e10; // invalid discriminant

    double sqrtD = sqrt(D);
    double z1 = (-b + sqrtD) / (2.0 * a);
    double z2 = (-b - sqrtD) / (2.0 * a);

    double z = NAN;
    if (z1 > 0.0 && R_finite(z1)) z = z1;
    else if (z2 > 0.0 && R_finite(z2)) z = z2;
    if (!R_finite(z) || z <= 0.0) return 1e10;

    double u = log(z);
    double Ae = alpha * z + (1.0 - alpha);
    double K = xprev * log(Ae) + lambda * (z - 1.0);
    double K2 = xprev * (alpha * (1.0 - alpha) * z) / (Ae * Ae) + lambda * z;
    if (K2 <= 0.0 || !R_finite(K2)) return 1e10;

    double logf_t = -0.5 * log(2.0 * M_PI * K2) + K - u * xt;
    if (!R_finite(logf_t)) return 1e10;

    ll += logf_t;
  }
  // return negative log-likelihood
  return -ll;
}
')

# ---------------------------------------------------------
# Wrapper R per ottimizzazione
# ---------------------------------------------------------
sp_inar1_poisson_rcpp <- function(x, init = NULL, control = list()) {
    if (!is.numeric(x) || length(x) < 2) stop("x must be numeric integer-like vector length >= 2")
    x <- as.integer(x)
    n <- length(x)

    if (is.null(init)) {
        mu <- mean(x)
        alpha0 <- 0.3
        lambda0 <- max(0.1, mu * (1 - alpha0))
        init <- c(alpha = alpha0, lambda = lambda0)
    }

    to_theta <- function(par) c(qlogis(par["alpha"]), log(par["lambda"]))
    from_theta <- function(theta) c(plogis(theta[1]), exp(theta[2]))

    nll_transformed <- function(theta) {
        par <- from_theta(theta)
        alpha <- par["alpha"]
        lambda <- par["lambda"]
        if (alpha <= 0 || alpha >= 1 || lambda <= 0) return(1e10)
        return(saddle_nll_cpp(x, alpha, lambda))
    }

    theta0 <- to_theta(init)
    opt <- optim(theta0, nll_transformed, method = "BFGS", control = control)
    theta_hat <- opt$par
    par_hat <- from_theta(theta_hat)

    list(par = par_hat, value = opt$value, conv = opt$convergence, message = opt$message)
}
# Parameter trick:
# 1. to_theta: par_transform: trasform parameters to unconstrained space
#    alpha in (0,1) -> logit(alpha) = log(alpha/(1-alpha)) in R
#    lambda in (0,inf) -> log(lambda) in R
# 2. from_theta: par_back: trasform back to original space
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



# Example
set.seed(42)
N = 1000
a1 = 0.3
l = 2
x <- genINAR(n, a1,l, inn = "poisson")$X
fit <- sp_inar1_poisson_rcpp(x)
print(fit$par)
print((fit$par - c(a1, l))^2)


