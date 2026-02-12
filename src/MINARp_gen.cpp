// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>

using namespace R;
using namespace Rcpp;

//' Generate a MINAR(p) series
//'
//' Experimental function
//'
//' @param resid Numeric matrix of size n x m, where n is the number of time points and m is the number of series.
//' @param A Numeric matrix of size m x m*p, where m is the number of series. Each column of A contains the autoregressive parameters for the corresponding lag.
//' @details
//' This is an internal function.
// [[Rcpp::export]]
NumericMatrix MINARp_gen_cpp(NumericMatrix resid, NumericMatrix A) {
    unsigned int n = resid.nrow();
    unsigned int m = resid.ncol();
    unsigned int p = A.ncol()/m;
    std::vector<NumericMatrix> Ai(p);
    int vals = 0;
    NumericMatrix sim(n,m);

    // divide A that is m x m*p into p squared matrices Ai of size m x m, each one containing the autoregressive parameters for the corresponding p lags
    for (unsigned int k = 0; k < p; k++) {
        Ai[k] = A(Range(0, m - 1), Range (k * m, (k + 1) * m - 1));
    }
    // print(Ai[0]);
    // print(Ai[1]);

    sim = clone(resid);
    for (unsigned int t = p; t < n; t++) {
        for (unsigned int j = 0; j < m; j++) {
            vals = 0;
            for (unsigned int k = 0; k < p; k++) {
                for (unsigned int h = 0; h < m; h++) {
                    vals += R::rbinom(sim(t - k - 1, h), Ai[k](h,j));
                }
            }
            sim(t, j) = vals + resid(t,j);
        }
    }

    return sim;
}

/*** R
# x <- rpois(1000,2.5)
# y <- rpois(1000,5)
# z <- rpois(1000,1)
# x <- c(0,3,4,5)
# y <- c(1,0,2,3)
# z <- c(0,5,6,0)
# RES <- cbind(x,y,z)
# A1 <- matrix(c(0.5,0.2,0.1,0.3,0.4,0.2,0.1,0.1,0.3),nrow=3,ncol=3)
# A2 <- matrix(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),nrow=3,ncol=3)
# A <- cbind(A1,A2)
# MINARp_gen_cpp(RES,A)
#
# set.seed(1234)
# x <- rpois(1000,2.5)
# y <- rpois(1000,5)
# z <- rpois(1000,1)
# RES <- cbind(x,y,z)
# A1 <- matrix(c(0.5,0,0,0,0.4,0,0,0,0.2),nrow=3,ncol=3)
# A2 <- matrix(c(0.2,0,0,0,0.3,0,0,0,0.4),nrow=3,ncol=3)
# A <- cbind(A1,A2)
# set.seed(1234)
# MINARp_gen_cpp(RES,A)
# set.seed(1234)
# MINARp_gen_cpp(RES,A) |> cor()
*/
