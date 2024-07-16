// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <PoissonBinomial.h>

using namespace R;
using namespace Rcpp;
using namespace PoissonBinomial;


//' Harris-McCabe score statistic to test for dependence in an integer autoregressive process
//' @param x NumericVector
//' @param method unsigned int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector HMC_Cpp(NumericVector x){
  int n = x.length();
  double mean_x = mean(noNA(x));
  double sd_x = sd(noNA(x));

  NumericVector out(2);
  Range idx = seq(1,n-1);
  Range idx_1 = seq(0,n-2);

  NumericVector  xsum = x[idx];
  NumericVector  xsum_1 = x[idx_1];

  IntegerVector tab = table(xsum);
  IntegerVector tab_1 = table(xsum - 1);

  NumericVector ftab = as<NumericVector>(tab) / (n - 1);
  NumericVector ftab_1 = as<NumericVector>(tab_1) / (n - 1);

  NumericVector pi_hat(n-1);
  NumericVector pi_hat_L1(n-1);
  for (int j = 0; j < n-1; ++j) {
      pi_hat[j] = ftab[xsum[j]];
      pi_hat_L1[j] = ftab_1[xsum[j]];

      //     std::cout << std::isnan(s_temp[i]) << std::endl;
      Rprintf("xsum[j] is: %f \n", xsum[j]);
      Rprintf("pi_hat[j] is: %f \n", pi_hat[j]);
      Rprintf("pi_hat_L1[j] is: %f \n", pi_hat_L1[j]);


      if(std::isnan(pi_hat[j])){
          pi_hat[j] = 0;
      }
      if(std::isnan(pi_hat_L1[j])){
          pi_hat_L1[j] = 0;
      }
  }


  NumericVector g_t = pi_hat_L1/pi_hat;
  // double mean_g = mean(g_t);
  double sd_g = sd(g_t);

  double NUM = sum((xsum_1-mean_x)*(g_t-1));

  double stat = NUM/(sd_x*sd_g);

  out[0] = stat/sqrt(n);
  out[1] = 1 - R::pnorm(out[0],0.0, 1.0, 1, 0);

  return out;
}



//' Semiparametric bootstrap version of the Harris McCabe score test.
//' @param x NumericVector
//' @param B int
//' @param method unsigned int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector HMC_semiparBOOT_Cpp(NumericVector x, int B){
    int n = x.length();
    unsigned int niter;

  NumericVector s_temp(B);
  for(int i = 0; i < B; i++){

    NumericVector xb(n);
    niter = 0;
    // print(xb);

    // check!
    LogicalVector id(n);
    do {
        id = xb==xb[0];

        while(Rcpp::all(id).is_true()) {
            xb = RcppArmadillo::sample(x,n,true);
            id = xb==xb[0];
        }

        s_temp[i] = HMC_Cpp(xb)[0];

        // asd = std::isnan(s_temp[i]);
        // if(asd){
        //     std::cout << std::isnan(s_temp[i]) << std::endl;
        //     // Rprintf("condition is: %d \n", asd);
        //     // Rprintf("s_temp is: %f \n", s_temp[i]);
        // }

        niter += 1;

        // con false esce, con true resta
    } while ( std::isnan(s_temp[i]) & (niter < 10) );

  }

  return s_temp;
}




// [[Rcpp::export]]
List HMCtest_boot(NumericVector X, unsigned int arrival, unsigned int type, int B){
    // int n = X.length();

    double Smc = HMC_Cpp(X)[0];
    double B_stat;
    double B_pval;
    NumericVector SmcB(B);

    type = 1;
    if(type == 1){
        // SEMIPARAMETRIC
        SmcB = HMC_semiparBOOT_Cpp(X,B);
    }
    // else if(type == 2){
    //     // PARAMETRIC
    //     SmcB = sunMC_parBOOT_Cpp(X,B,arrival);
    // }
    // else if(type == 3){
    //     // PIT
    //     SmcB = sunMC_pitBOOT_Cpp(X,B,arrival);
    // }

    B_stat = mean(SmcB);
    B_pval = mean(abs(SmcB) > std::fabs(Smc));


    return List::create(
        _["stat"]  = B_stat,
        _["pval"]  = B_pval
    );
}

/*** R
# x <- rpois(5000,2)
# B <- 1000
# method <- 1
# # test 1:
# S <- HMC_Cpp(x,method)
# Sb <- HMC_semiparBOOT_Cpp(x,B,method)
# mean(S[1] > Sb) # pval boot
# S[2] # pval normale
# # test 2: FUNZIONE FINALE WRAPPED
# arrival <- 1
# HMCtest_boot(x,arrival,method,B)
*/
