// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

using namespace R;
using namespace Rcpp;

// source :
// https://teuder.github.io/rcpp4everyone_en/220_dpqr_functions.html


// // [[Rcpp::export]]
// double TAsum_cpp(double sim, double arrivals, double a) {
//   double TAsum;
//   double thinning;
//
//   thinning = R::rbinom(sim,a);
//   TAsum = thinning + arrivals;
//
//   return TAsum;
// }


// [[Rcpp::export]]
NumericVector INAR1_cpp(NumericVector resid, double a) {

  unsigned int n = resid.length();
  NumericVector sim(n);

  sim[0] = resid[0];
  for (unsigned int i = 1; i < n; i++) {
    sim[i] = R::rbinom(sim[i-1],a) + resid[i];
  }

  return sim;
}

/*** R
# aa <- 0.8
# x <- rpois(1000,2)
# y <- INAR1_cpp(x,aa)
# par(mfrow=c(2,1))
# plot(x,type="b")
# plot(y,type="b")
# par(mfrow=c(1,1))
# acf(y)$acf[2]
*/
