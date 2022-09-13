#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//


// a_ <- unname(par[1]) # alpha
//   for(i in 2:s){
//     sim_X[i] <- rbinom(1,sim_X[i-1],prob=a_) + resid_[i]
//   }
//   dati_sim <- data.frame(X=sim_X[(burnout+1):s],res=resid_[(burnout+1):s])


// [[Rcpp::export]]
IntegerVector rcpp_genINARmarg(double a,IntegerVector Xseries,IntegerVector resid) {

  // n=lunghezza del vettore resid
  int n = resid.size();

  // genero vettore di lunghezza n e con tutti 0
  // sti -9999 saranno tuttti sovrascritti, tranne il primo!
  IntegerVector out(n,-9999);

  // e mo se chiagne
  // For loop, nota che cpp index shift to 0 !!!PER VETTORI!!!
  for (int i=1; i <= n; ++i) {
  //   // rbinom(n, size, prob) è NumericVector
  //   NumericVector rbin = Rcpp::rbinom(1, 5, a);
  int rbin = as<int>(rbinom(1, Xseries[i-1], a));
    out[i] = rbin + resid[i];
  }
  // prendo solo valori non-negativi, tolgo il primo el. inutile
    return out[out >=0];
}
