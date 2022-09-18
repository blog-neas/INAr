// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

using namespace R;
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector INAR1_cpp(NumericVector resid, double a) {

    unsigned int n = resid.length();
    NumericVector sim(n);

    sim = clone(resid);
    for (unsigned int i = 1; i < n; i++) {
        sim[i] = R::rbinom(sim[i-1],a) + resid[i];
    }

    return sim;
}

// [[Rcpp::export]]
NumericVector INARp_cpp(NumericVector resid, DoubleVector a) {

    unsigned int n = resid.length();
    unsigned int lags = a.length();
    int vals = 0;
    NumericVector sim(n);

    sim = clone(resid);
    for (unsigned int i = lags; i < n; i++) {
        vals = 0;
        for(unsigned int j = 0 ; j < lags; j++) {
            vals += R::rbinom(sim[i - lags + j],a[j]);
        }
        sim[i] = vals + resid[i];
        // sim[i] = R::rbinom(sim[i-1],a) + resid[i];
    }

    return sim;
}

/*** R
# # check INAR(p)
# aa <- c(0.1,0.3,0.1,0.2,0.2)
# x <- rpois(1000,2)
# y <- INARp_cpp(x,aa)
# par(mfrow=c(2,1))
# plot(x,type="b")
# plot(y,type="b")
# par(mfrow=c(1,1))
# acf(y,plot = FALSE)$acf[2]
# # check INAR(1)
# aa <- 0.8
# x <- rpois(1000,2)
# y <- INAR1_cpp(x,aa)
# par(mfrow=c(2,1))
# plot(x,type="b")
# plot(y,type="b")
# par(mfrow=c(1,1))
# acf(y,plot = FALSE)$acf[2]
# # check if INAR1_cpp and INARp_cpp generate the same series
# aa <- 0.8
# x <- rpois(1000,2)
# set.seed(1234)
# y1 <- INAR1_cpp(x,aa)
# set.seed(1234)
# y2 <- INARp_cpp(x,aa)
# plot(y1,type="b",col="blue")
# lines(y2,type="l",col="red")
*/
