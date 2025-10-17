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
