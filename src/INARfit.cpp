#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace R;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat YW_cpp(const arma::vec& r) {
    int p = r.n_elem;
    arma::mat R = arma::eye(p, p);

    for(int i = 0; i < p - 1; i++) {
        for(int j = i + 1; j < p; j++) {
            // Adjust index: subtract 1 to convert to 0-indexing
            R(i, j) = r(j - i - 1);
        }
    }
    // upper symmetrize
    R = arma::symmatu(R);
    return R;
}


// [[Rcpp::export]]
NumericVector INARfitted_cpp(NumericVector X, double mINN, DoubleVector a) {

    unsigned int n = X.length();
    unsigned int p = a.length();
    double vals = 0.0;
    NumericVector fitted = clone(X);

    // calcolo i fitted senza sovrascrivere i lag osservati di X
    for (unsigned int t = p; t < n; t++) {
        vals = 0.0;
        for(unsigned int k = 0 ; k < p; k++) {
            vals += X[t - k - 1] * a[k];
        }
        fitted[t] = vals + mINN;
    }

    return fitted;
}

//[[Rcpp::export]]
Rcpp::List INARforecast_cpp(NumericVector X, double mINN, DoubleVector a, int h, int B, double alpha = 0.95) {

    unsigned int n = X.length();
    unsigned int p = a.length();
    double vals = 0.0;
    double eps = 0.0;
    NumericVector forecast(h);
    NumericVector lower(h);
    NumericVector upper(h);
    NumericMatrix forecastBoot(h,B);

//     calcolo i forecast
    if(B==0){
        for (unsigned int t = n; t < n + h; t++) {
            vals = 0.0;
            for(unsigned int k = 0 ; k < p; k++) {
                if(t - k - 1 < n) {
                    // CHECK come fare forecast giusto (h-step ahead))
                    // vals += R::dbinom(X[t - k - 1], 1, a[k]) * a[k]; // E[Binomial(X[t-k-1], a[k])]
                    vals += X[t - k - 1] * a[k];
                } else {
                    // CHECK come fare forecast giusto (h-step ahead)
                    // vals += R::dbinom(forecast[t - n - k - 1], 1, a[k]) * a[k]; // E[Binomial(forecast[t-k-1], a[k])]
                    vals += forecast[t - n - k - 1] * a[k];
                }
            }
            forecast[t - n] = vals + mINN;
        }
    } else {
        // bootstrap forecast
        // IDEA SCEMA, VEDI BISAGLIA GEROLIMETTO PER APPROFONDIMENTI:
        // faccio B simulazioni della serie con INARp_cpp, e per ogni sim
        // calcolo i forecast h-step ahead,
        // poi faccio la media dei forecast ottenuti dalle B simulazioni
        for (int b = 0; b < B; b++) {
            // NumericVector sim = INARp_cpp(X, a);
            for (unsigned int t = n; t < n + h; t++) {
                vals = 0.0;
                for(unsigned int k = 0 ; k < p; k++) {
                    if(t - k - 1 < n) {
                        // vals += sim[t - k - 1] * a[k];
                        vals += R::rbinom(X[t - k - 1], a[k]);
                    } else {
                        // vals += forecast[t - n - k - 1] * a[k];
                        vals += R::rbinom(forecastBoot(t - n - k - 1,b), a[k]);
                    }
                }

                // SELETTOPRE PER GENERARE INN DA DISTRIBUZIONE SCELTA
                // TO DO: implementare anche il caso semiparametrico, che è più complesso
                // if(inncode == 1){
                //     // Poisson case
                //     mINN = R::rpois(mINN);
                // }else if(inncode == 2){
                //     // NegBin case
                //     mINN = R::rnbinom(mINN, 1.0/(1.0+mINN));
                // }else if(inncode == 3){
                //     // GenPoi case
                //     // TO DO!
                //     mINN = R::rpois(mINN);
                // }else if(inncode == 4){
                //     // semiparametric case
                //     // TO DO!
                //     mINN = ;
                // }

                // per adesso solo POI-INAR(1)
                eps = R::rpois(mINN);
                forecastBoot(t - n, b) = vals + eps;
            }
        }
        forecast = rowMeans(forecastBoot);
    }
    return List::create(
        _["forecast"] = forecast,
        _["paths"] = forecastBoot
    );
}

