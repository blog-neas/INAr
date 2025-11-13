#include <Rcpp.h>
using namespace Rcpp;

// Poisson INAR(1) log-likelihood computation for ML estimation
// [[Rcpp::export]]
double inar1_poi_loglik_cpp(const IntegerVector& x, double alpha, double lambda) {
    int n = x.size();
    if (n < 2) return NA_REAL;
    double ll = 0.0;

    for (int t = 1; t < n; t++) {
        int xt = x[t];
        int xprev = x[t - 1];
        if (xt < 0 || xprev < 0) return NA_REAL;

        int kmax = std::min(xt, xprev);
        std::vector<double> log_p(kmax + 1);

        for (int k = 0; k <= kmax; k++) {
            // combinazione log
            double lchoose_val = R::lchoose(xprev, k);
            // log densità Poisson (x_t - k)
            int diff = xt - k;
            if (diff < 0) { log_p[k] = -INFINITY; continue; }

            double lp = lchoose_val + k * log(alpha) +
                (xprev - k) * log(1.0 - alpha) +
                R::dpois(diff, lambda, true);
            log_p[k] = lp;
        }

        // log-sum-exp trick per stabilità numerica
        double m = -INFINITY;
        for (int k = 0; k <= kmax; k++) if (log_p[k] > m) m = log_p[k];
        double sum_exp = 0.0;
        for (int k = 0; k <= kmax; k++) sum_exp += exp(log_p[k] - m);
        double log_p_xt = m + log(sum_exp);

        if (!R_finite(log_p_xt)) return -INFINITY;
        ll += log_p_xt;
    }
    return ll;
}

// Poisson INAR(1) log-likelihood computation for SP estimation
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
    // return positive log-likelihood (put -ll on the outside for minimization)
    return ll;
}

