// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
// #include <PoissonBinomial.h>
// #include "INAr.h"

using namespace R;
using namespace Rcpp;
// using namespace RcppParallel;
// using namespace PoissonBinomial;

//' Autocorrelation rho test statistic to test for unit root in an integer autoregressive process
//' @param x NumericVector
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector RHO_Cpp(NumericVector x){
     int n = x.length();
     NumericVector out(2);
     int lag = 1; // to compute acf at lag 1

     NumericVector x_scaled = x - mean(noNA(x));
     NumericVector xsum = x_scaled[Range(lag,n-1)];
     NumericVector xsum_1 = x_scaled[Range(0,n-lag-1)];

     double NUM = sum(xsum*xsum_1);
     double DEN = sum(pow(x_scaled,2));

     double stat = NUM/DEN; // rho_1

     out[0] = stat * sqrt(n);
     out[1] = 1 - R::pnorm(out[0],0.0, 1.0, 1, 0);

     return out;
}


//' Autocorrelation bootstrap test.
//' @param x NumericVector
//' @param B int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector RHO_BOOT_Cpp(NumericVector x, int B){
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

             s_temp[i] = RHO_Cpp(xb)[0];

             niter += 1;

             // con false esce, con true resta
         } while ( std::isnan(s_temp[i]) & (niter < 10) );

     }

     return s_temp;
}



//
//  ----------------------------------------------------------------------------
//  ----------------------------------------------------------------------------
//  ----------------------------------------------------------------------------
//

// testing the possibility to use a parallel version of the bootstrap
// WARNING: compilation error on Mac OS X
//
// // [[Rcpp::plugins(openmp)]]
// #include <omp.h>
//
// // [[Rcpp::export]]
// NumericVector RHO_BOOT_Cpp_Parallel(NumericVector x, int B, int num_threads = 2) {
//     int n = x.length();
//     NumericVector s_temp(B);
//
//         // Set the number of threads (optional, default is usually the number of cores)
//         if(num_threads > 1 & num_threads <= omp_get_max_threads()){
//             omp_set_num_threads(num_threads);
//         }else{
//             num_threads = 1;
//         }
//
// #pragma omp parallel num_threads(num_threads)
// {
//     // Each thread needs its own RNG
//     RNGScope scope;
//
// #pragma omp for
//     for(int i = 0; i < B; i++) {
//         // simplified version for testing, to update using the correct do while loop
//         NumericVector xb = RcppArmadillo::sample(x, n, true);
//         s_temp[i] = RHO_Cpp(xb)[0];
//     }
// }
//
// return s_temp;
// }


// #include <omp.h>
//
// // Placeholder for the RHO_Cpp function, which you need to provide
// NumericVector RHO_Cpp(NumericVector x);
//
// // [[Rcpp::export]]
// NumericVector RHO_BOOT_Cpp_Parallel(NumericVector x, int B, int threads = 0) {
//     int n = x.length();
//     NumericVector s_temp(B);
//
//     // Set the number of threads (optional, default is usually the number of cores)
//     if(threads > 0 & threads <= omp_get_max_threads()){
//         omp_set_num_threads(threads);
//     }else{
//         omp_set_num_threads(1);
//     }
//
// #pragma omp parallel for
//     for(int i = 0; i < B; i++) {
//         NumericVector xb = RcppArmadillo::sample(x, n, true);
//         s_temp[i] = RHO_Cpp(xb)[0];
//     }
//
//     return s_temp;
// }


//
//  ----------------------------------------------------------------------------
//  ----------------------------------------------------------------------------
//  ----------------------------------------------------------------------------
//


//
// PARTE Script Sun-McCabe, prima era in INARbootSMC.cpp -----
//

//' Sun-McCabe score statistic to test for dependence in an integer autoregressive process
//' @param x NumericVector
//' @param method unsigned int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector SMC_Cpp(NumericVector x, unsigned int method){
  int n = x.length();

  NumericVector out(2);
  Range idx = seq(1,n-1); // delete and use Range(..)
  Range idx_1 = seq(0,n-2); // delete and use Range(..)

  NumericVector xsum = x[idx]; // change with x[Range(1,n-1)]
  NumericVector xsum_1 = x[idx_1]; // change with x[Range(0,n-2)]

  if(method == 1){
    // poisson case

    double mu_x = mean(noNA(x)); // = lambda

    NumericVector xsumL = xsum - mu_x;
    NumericVector xsumL1 = xsum_1 - mu_x;
    NumericVector NUM = xsumL*xsumL1;

    double stat = sum(NUM)/mu_x;

    out[0] = stat/sqrt(n);
    out[1] = 1 - R::pnorm(out[0],0.0, 1.0, 1, 0);
  }

  if(method == 2){
    // negbin case

    double mu_x = mean(noNA(x));
    double var_x = var(noNA(x));
    double sd_x = sd(noNA(x));
    double diffvarmu = std::fabs(var_x - mu_x); // # trick, uso VAL ASS DIFF

    // check underdispersion
    if(var_x <= mu_x){
        // Rcpp::Rcout << "Underdispersion detected: switching to generalized Poisson method.\n";
        // method = 3;
        out[0] = NAN;
        out[1] = NAN;
        return out;
    }
    // printf("1) %f, %f %f \n",mu_x,var_x,std::abs(var_x-mu_x));
    // printf("2) %f, %f \n",pow(mu_x,2),var_x-mu_x);

    // CHECK!
    double r_hat = pow(mu_x,2)/diffvarmu;
    double p_hat = diffvarmu/var_x;

    NumericVector ind_val = ifelse(xsum > 0, xsum/((r_hat+xsum-1)*p_hat), 0);

    // for(int i=0; i<asd.length(); ++i){
    //   Rprintf("the value of v[%i] : %f | xsum[%i] : %f  | (%f) \n", i, asd[i], i, xsum[i],xsum[i]/((r_hat+xsum[i]-1)*p_hat));
    // }

    NumericVector xsumL = ind_val - 1;
    NumericVector xsumL1 = xsum_1 - mu_x;
    NumericVector NUM = xsumL*xsumL1;

    double sd_g = sd(xsumL); // era sd(ghat)

    double stat = sum(NUM)/(sd_x*sd_g);
    out[0] = stat/sqrt(n);
    out[1] = 1 - R::pnorm(out[0],0.0, 1.0, 1, 0);
  }
  if(method == 3){
    // genpoi case

    double mu_x = mean(noNA(x));
    // double var_x = var(noNA(x));
    double sd_x = sd(noNA(x));

    double lambda_hat = sqrt(pow(mu_x,3))/sd_x; // # trick
    double kappa_hat = 1 - sqrt(mu_x)/sd_x; // # trick
    // printf("%f, %f \n",lambda_hat,kappa_hat);

    // NumericVector xsumL = ind_val - 1;
    NumericVector xsumL1 = xsum_1 - mu_x;

    NumericVector ltmp1 = (xsum - 2)*log(lambda_hat + kappa_hat*(xsum-1));
    NumericVector ltmp2 = (xsum - 1)*log(lambda_hat + kappa_hat*(xsum));
    NumericVector tmp3 = exp(kappa_hat)*xsum;

    NumericVector g_hat = exp(ltmp1-ltmp2)*tmp3 - 1;
    double mu_g = mean(g_hat);
    double sd_g = sd(g_hat);

    // for(int i=0; i < ltmp1.length(); ++i){
    //   Rprintf("the value of tmp1[%i] : %f \n", i, ltmp1[i]);
    //   Rprintf("the value of tmp2[%i] : %f \n", i, ltmp2[i]);
    //   Rprintf("the value of tmp3[%i] : %f \n", i, tmp3[i]);
    //   Rprintf("the value of ghat[%i] : %f \n", i, g_hat[i]);
    // }

    NumericVector xsumL = g_hat - mu_g;

    NumericVector NUM = xsumL*xsumL1;

    double stat = sum(NUM)/(sd_x*sd_g);
    out[0] = stat/sqrt(n);
    out[1] = 1 - R::pnorm(out[0],0.0, 1.0, 1, 0);
  }
  if(method == 4){
    // katz case
    // in development

    double mu_x = mean(noNA(x));
    double var_x = var(noNA(x));
    double sd_x = sd(noNA(x));

    double a_hat = pow(mu_x,2)/var_x;
    double b_hat = 1 - mu_x/var_x;

    NumericVector ind_val = ifelse(xsum > 0, xsum/(a_hat-b_hat+b_hat*xsum), 0);

    NumericVector g_hat = ind_val - 1;
    double mu_g = mean(g_hat);
    double sd_g = sd(g_hat);

    NumericVector xsumL = g_hat - mu_g;
    NumericVector xsumL1 = xsum_1 - mu_x;
    NumericVector NUM = xsumL*xsumL1;

    double stat = sum(NUM)/(sd_x*sd_g);
    out[0] = stat/sqrt(n);
    out[1] = 1 - R::pnorm(out[0],0.0, 1.0, 1, 0);
  }
  return out;
}



//' Semiparametric bootstrap version of the Sun-McCabe score test.
//' @param x NumericVector
//' @param B int
//' @param method unsigned int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector SMC_semiparBOOT_Cpp(NumericVector x, int B, unsigned int method){
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

            s_temp[i] = SMC_Cpp(xb,method)[0];

            niter += 1;

            // con false esce, con true resta
        } while ( std::isnan(s_temp[i]) & (niter < 10) );

    }

    return s_temp;
}


//' Parametric bootstrap version of the Sun-McCabe score test.
//' @param x NumericVector
//' @param B int
//' @param method unsigned int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector SMC_parBOOT_Cpp(NumericVector x, int B, unsigned int method){
    int n = x.length();
    unsigned int niter;
    double mean_x = mean(noNA(x));
    double var_x = var(noNA(x));

    // Function genUniGpois("GenUniGpois");

    // NumericVector out(3);
    // NumericVector MB(B);
    // NumericVector VB(B);

    NumericVector s_temp(B);
    if(method==1){

        // POISSON
        for(int i = 0; i < B; i++){
            niter = 0;
            NumericVector xb(n);
            LogicalVector id(n);
            do {
                id = xb==xb[0];

                while(Rcpp::all(id).is_true()) {
                    // SOLO CASO POISSON
                    // lambda_x = mean_x
                    xb = Rcpp::rpois(n,mean_x);
                    id = xb==xb[0];
                }

                s_temp[i] = SMC_Cpp(xb,method)[0];

                niter += 1;
                // con false esce, con true resta
            } while ( std::isnan(s_temp[i]) & (niter < 10) );
        }
    }
    if(method==2){
        // NEGBIN
        for(int i = 0; i < B; i++){
            NumericVector xb(n);
            // print(xb);

            // check!
            LogicalVector id(n);
            id = xb==xb[0];
            niter = 0;

            do {
                id = xb==xb[0];

                while(Rcpp::all(id).is_true()) {
                    // SOLO CASO NEGBIN
                    // secondo formula (22) Sun-McCabe, Katz applicata a NegBin:
                    // gamma_x = mean_x^2/(var_x-mean_x);
                    // p_x = 1 - mean_x/var_x = (var_x - mean_x)/var_x;
                    // p_x = abs(var_x - mean_x)/var_x ; // TRICK!
                    // pcompl_x = 1-p_x = mean_x/var_x;
                    double diffvarmu = std::fabs(var_x - mean_x); // # trick, uso VAL ASS DIFF
                    double gamma_x = pow(mean_x,2)/diffvarmu;
                    double pcompl_x = 1 - diffvarmu/var_x;
                    xb = Rcpp::rnbinom(n, gamma_x, pcompl_x);
                    // // controllo underdispersion
                    // if(var(xb) <= mean(xb)){
                    //     xb = rep(0,n);
                    // }
                    id = xb==xb[0];
                }

                // controllo implementato in suMC_Cpp
                // if(var(xb) <= mean(xb)){
                // s_temp[i] = nan();
                // }
                s_temp[i] = SMC_Cpp(xb,method)[0];


                niter += 1;
                // con false esce, con true resta
            } while ( std::isnan(s_temp[i]) & (niter < 10) );
        }
    }
    if(method==3){
        // GENPOI
        for(int i = 0; i < B; i++){
            niter = 0;
            NumericVector xb(n);
            LogicalVector id(n);
            do {
                id = xb==xb[0];


                while(Rcpp::all(id).is_true()) {
                    //
                    double lambda = sqrt(pow(mean_x,3)/var_x); // is theta in genUniGpois()
                    double kappa = 1 - sqrt(mean_x/var_x); // is lambda in genUniGpois()
                    double w = exp(-kappa);
                    // NumericVector myset(n);
                    for(int j = 0; j < n; j++){
                        double mys = exp(-lambda);
                        double myp = mys;
                        int xx = 0;
                        double u = R::runif(0,1);
                        while(u > mys){
                            xx++;
                            double myc = lambda + kappa * (xx - 1);
                            myp = w * myc * pow(1 + kappa/myc, xx - 1) * myp * pow(xx,-1);
                            mys = mys + myp;
                        }
                        xb[j] = xx;
                    }
                    id = xb==xb[0];
                }

                s_temp[i] = SMC_Cpp(xb,method)[0];

                niter += 1;
                // con false esce, con true resta
            } while ( std::isnan(s_temp[i]) & (niter < 10) );
        }

    }

    return s_temp;
}


// //' Function to generate generalized Poisson random variables
//  //' @param n int
//  //' @param lambda double
//  //' @param kappa double
//  //' @details
//  //' This is an internal function, it will be excluded in future versions.
//  //' @export
//  // [[Rcpp::export]]
// NumericVector rgenpois(int n, double lambda, double kappa){
//     NumericVector result(n);
//     double w = exp(-kappa);
//     for(int j = 0; j < n; j++){
//         double mys = exp(-lambda);
//         double myp = mys;
//         int xx = 0;
//         double u = R::runif(0,1);
//         while(u > mys){
//             xx++;
//             double myc = lambda + kappa * (xx - 1);
//             myp = w * myc * pow(1 + kappa/myc, xx - 1) * myp * pow(xx,-1);
//             mys = mys + myp;
//         }
//         result[j] = xx;
//     }
//     return result;
// }



//' PIT bootstrap version of the Sun-McCabe score test.
//' @param x NumericVector
//' @param B int
//' @param method unsigned int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' !!!!! DA IMPLEMENTARE !!!!!
//' @export
// [[Rcpp::export]]
NumericVector SMC_pitBOOT_Cpp(NumericVector x, int B, unsigned int method){
    int n = x.length();
    unsigned int niter;
    double mean_x = mean(noNA(x));
    double var_x = var(noNA(x));

    NumericVector s_temp(B);
    if(method==1){
        // POISSON
        for(int i = 0; i < B; i++){
            NumericVector ub(n);
            NumericVector xb(n);
            DoubleVector q(n);
            DoubleVector u(n);
            LogicalVector id(n);

            niter = 0;
            do {
                id = xb==xb[0];

                while(Rcpp::all(id).is_true()) {
                    // SOLO CASO POISSON
                    // lambda_x = mean_x
                    // xb = Rcpp::rpois(n,mean_x);
                    // id = xb==xb[0];

                    q = Rcpp::runif(n);
                    u = q*Rcpp::ppois(x, mean_x) + (1-q)*Rcpp::ppois(x-1, mean_x);

                    ub = RcppArmadillo::sample(u,n,true);
                    xb = Rcpp::qpois(ub, mean_x);
                    id = xb==xb[0];

                }

                s_temp[i] = SMC_Cpp(xb,method)[0];

                niter += 1;
                // con false esce, con true resta
            } while ( std::isnan(s_temp[i]) & (niter < 10) );
        }
    }
    if(method==2){
        // NEGBIN
        for(int i = 0; i < B; i++){
            NumericVector ub(n);
            NumericVector xb(n);
            DoubleVector q(n);
            DoubleVector u(n);
            LogicalVector id(n);

            niter = 0;
            do {
                id = xb==xb[0];

                while(Rcpp::all(id).is_true()) {
                    // SOLO CASO NEGBIN
                    double diffvarmu = std::fabs(var_x - mean_x); // # trick, uso VAL ASS DIFF
                    double gamma_x = pow(mean_x,2)/diffvarmu;
                    double pcompl_x = 1 - diffvarmu/var_x;

                    q = Rcpp::runif(n);
                    u = q*Rcpp::pnbinom(x, gamma_x, pcompl_x) + (1-q)*Rcpp::pnbinom(x-1, gamma_x, pcompl_x);

                    ub = RcppArmadillo::sample(u,n,true);
                    xb = Rcpp::qnbinom(ub, gamma_x, pcompl_x);
                    id = xb==xb[0];
                }

                // controllo implementato in suMC_Cpp
                // if(var(xb) <= mean(xb)){
                // s_temp[i] = nan();
                // }
                s_temp[i] = SMC_Cpp(xb,method)[0];


                niter += 1;
                // con false esce, con true resta
            } while ( std::isnan(s_temp[i]) & (niter < 10) );

        }
    }

    return s_temp;
}



//' Wrapper function for compution the Sun-McCabe bootstrap score test.
//' @param X NumericVector
//' @param arrival int
//' @param type unsigned int
//' @param B int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
List SMCtest_boot(NumericVector X, unsigned int arrival, unsigned int type, int B){
    // int n = X.length();

    double Smc = SMC_Cpp(X,arrival)[0];
    double B_stat;
    double B_pval;
    NumericVector SmcB(B);

    if(type == 1){
        // SEMIPARAMETRIC
        SmcB = SMC_semiparBOOT_Cpp(X,B,arrival);
    }
    else if(type == 2){
        // PARAMETRIC
        SmcB = SMC_parBOOT_Cpp(X,B,arrival);
    }
    else if(type == 3){
        // PIT
        SmcB = SMC_pitBOOT_Cpp(X,B,arrival);
    }

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
# S <- SMC_Cpp(x,method)
# Sb <- SMC_semiparBOOT_Cpp(x,B,method)
# mean(S[1] > Sb) # pval boot
# S[2] # pval normale
# # test 2: FUNZIONE FINALE WRAPPED
# arrival <- 1
# SMCtest_boot(x,arrival,method,B)
#
# x <- rpois(100,3)
# Pb <- SMC_parBOOT_Cpp(x,99,3)
*/



//
// PARTE Script Harris-McCabe, prima era in INARbootHMC.cpp -----
//



//' Sort and remove duplicates from a numeric vector
//' @param v NumericVector
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
//[[Rcpp::export]]
NumericVector sortunique(NumericVector v) {
     NumericVector sv = Rcpp::unique(v);
     std::sort(sv.begin(), sv.end());
     return sv;
}

//' Compute the empirical cumulative distribution function of a numeric vector
//' @param eval NumericVector
//' @param x NumericVector
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
//[[Rcpp::export]]
DataFrame ecdfcpp(NumericVector eval, NumericVector x) {
     // Valori unici ordinati in eval
     NumericVector samp = sortunique(eval);

     int n = samp.size();
     int m = x.size();
     int count = 0;

     // Vettore per le frequenze cumulate
     IntegerVector abs(n);
     // Vettore per le probabilità puntuali
     NumericVector probs(n);
     // Vettore per le frequenze cumulate relative
     NumericVector relcum(n);
     // Vettore per la probabilità di fallback
     NumericVector fallback_probs(n);

     // Calcola la probabilità relativa cumulata e la probabilità di fallback
     double min_prob = probs[0];
     relcum[0] = probs[0];

     for (int i = 0; i < n; ++i){
         // frequenze assolute cumulate step
         abs[i] = std::lower_bound(x.begin(), x.end(), samp[i]) - x.begin();
         count = 0;
         for (int j = 0; j < m; ++j){
             if (x[j] == samp[i]) {
                 count++;
             }
         }

         // frequenze relative cumulate step
         probs[i] = (double)count / m;

         // relcum step
         relcum[i] = relcum[i - 1] + probs[i];
         // fallback step
         if (probs[i] == 0) {
             fallback_probs[i] = min_prob;
         } else {
             fallback_probs[i] = probs[i];
         }
         if (probs[i] > 0) {
             min_prob = probs[i];
         }
     }

     return DataFrame::create(
         Named("value") = samp,
         Named("abscum") = abs,
         Named("relcum") = relcum,
         Named("probability") = probs,
         Named("fallback") = fallback_probs
     );
}


//' Harris-McCabe score statistic to test for dependence in an integer autoregressive process
//' @param x NumericVector
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector HMC_Cpp(NumericVector x){
     int n = x.length();
     double mean_x = mean(noNA(x));
     double sd_x = sd(noNA(x));
     LogicalVector id(n);
     LogicalVector id_1(n);

     NumericVector out(2);
     Range idx = seq(1,n-1);
     Range idx_1 = seq(0,n-2);

     NumericVector  xsum = x[idx];
     NumericVector  xsum_1 = x[idx_1];

     // NumericVector x2 = xsum-1;
     NumericVector x2 = x-1;

     NumericVector x_full(x.size() + x2.size());

     // Merge the vectors using the merge function
     std::merge(x.begin(), x.end(), x2.begin(), x2.end(),
                x_full.begin());
     // std::cout << x_full << std::endl;

     NumericVector eval = sortunique(x_full);
     // std::cout << eval << std::endl;
     DataFrame TAB = ecdfcpp(eval, x);
     NumericVector values = TAB["value"];
     NumericVector relfreq = TAB["fallback"]; // improved estimatof probs. (check!)

     NumericVector pi_hat(n-1);
     NumericVector pi_hat_L1(n-1);

     for (int j = 0; j < n-1; j++){

         id = values==xsum[j];
         id_1 = values==xsum[j]-1;
         // std::cout << id << std::endl;

         NumericVector tmp = relfreq[id];
         // std::cout << tmp << std::endl;

         pi_hat[j] = as<double>(relfreq[id]);
         pi_hat_L1[j] = as<double>(relfreq[id_1]);
         // Rprintf("x[j] is: %f, pi_hat[j] is: %f, pi_hat_L1[j] is: %f and the ratio is %f \n", x[j], pi_hat[j], pi_hat_L1[j], pi_hat_L1[j]/pi_hat[j]);
     }

     NumericVector g_t = pi_hat_L1/pi_hat;
     // std::cout << x << std::endl;
     // std::cout << pi_hat_L1 << std::endl;
     // std::cout << pi_hat << std::endl;
     // std::cout << g_t << std::endl;

     double sd_g = sd(g_t);

     NumericVector xsumL = xsum_1 - mean_x;

     double NUM = sum(xsumL*(g_t-1));
     double stat = NUM/(sd_x*sd_g);

     out[0] = stat/sqrt(n);
     out[1] = 1 - R::pnorm(out[0],0.0, 1.0, 1, 0);
     // std::cout << out[0] << std::endl;

     return out;
}


//' Bootstrap version of the Harris-McCabe score test.
//' !!!WARNING!!! Still under development, do not use! It will be replaced by INARtest() in future versions.
//' @param x NumericVector
//' @param B int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
NumericVector HMC_BOOT_Cpp(NumericVector x, int B){
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

             if(!arma::is_finite(s_temp[i])){
                 s_temp[i] = sign(s_temp[i])*pow(n,10);
             };

             // asd = std::isnan(s_temp[i]);
             // if(asd){
             //     std::cout << std::isnan(s_temp[i]) << std::endl;
             //     // Rprintf("condition is: %d \n", asd);
             //     // Rprintf("s_temp is: %f \n", s_temp[i]);
             // }

             niter += 1;
             // con false esce, con true resta
             // check for whether a value is finite, e.g. not NaN,Inf, or -Inf, by using arma::is_finite()
         } while (!arma::is_finite(s_temp[i]) & (niter < 10) );

     }
     // std::cout << s_temp << std::endl;
     return s_temp;
}

//' Wrapper function for compution the Harris-McCabe bootstrap score test.
//' !!!WARNING!!! Still under development, do not use! It will be replaced by INARtest() in future versions.
//' @param X NumericVector
//' @param B int
//' @details
//' This is an internal function, it will be excluded in future versions.
//' @export
// [[Rcpp::export]]
List HMCtest_boot(NumericVector X, int B){
    // old input: unsigned int type
    // int n = X.length();

    double Smc = HMC_Cpp(X)[0];
    double B_stat;
    double B_pval;
    NumericVector SmcB(B);

    unsigned int type = 1;
    if(type == 1){
        SmcB = HMC_BOOT_Cpp(X,B);
    }
    // else if(type == 2){
    //     // PARAMETRIC
    //     SmcB = HMC_parBOOT_Cpp(X,B,arrival);
    // }
    // else if(type == 3){
    //     // PIT
    //     SmcB = HMC_pitBOOT_Cpp(X,B,arrival);
    // }

    B_stat = mean(SmcB);
    B_pval = mean(abs(SmcB) > std::fabs(Smc));


    return List::create(
        _["stat"]  = B_stat,
        _["pval"]  = B_pval
    );
}


/*** R
# set.seed(1913)
# x <- rpois(500,1)
# HMC_Cpp(x)
# x <- c(1, 2, 12, 0, 2, 3, 2, 3, 4, 0, 1, 1, 2, 2, 4)
# HMC_Cpp(seq(0,20,by=2))
# B <- 21
# # test 1:
# S <- HMC_Cpp(x)
# set.seed(1432)
# Sb <- HMC_semiparBOOT_Cpp(x,B)
# mean(S[1] > Sb) # pval boot
# S[2] # pval normale
# # test 2: FUNZIONE FINALE WRAPPED
# set.seed(1432)
# HMCtest_boot(x,B)
# ecdf_cpp(1,x);x
# cpplb(0,x);x
# ecdf_cpp(-3,x);sum(x <= -3)
# ecdf_cpp(5,x)/length(x);sum(x <= 5)/length(x)
# ecdf_cpp(5,x) - ecdf_cpp(4,x)
# cpplb(unique(sort(seq(-1,20))),x)
# ecdfcpp(seq(-1,16),x)
# ecdfcpp(unique(x),x)
#
# x <- c(1.1, 2.2, 3.3, 2.2, 1.1, 4.4, 3.3, 2.2)
# eval <- c(1.1, 2.2, 3.3, 5.5)
#
# # Utilizzo della funzione ecdfcpp
# ecdfcpp(eval, x)
# x <- rpois(100,2)
# microbenchmark::microbenchmark(
#     RHO_BOOT_Cpp(x,399),
#     RHO_BOOT_Cpp_Parallel(x,399),
#     RHO_BOOT_Cpp_Parallel(x,399,2),
#     RHO_BOOT_Cpp_Parallel(x,399,3),
#     RHO_BOOT_Cpp_Parallel(x,399,4),
#     RHO_BOOT_Cpp_Parallel(x,399,5),
#     RHO_BOOT_Cpp_Parallel(x,399,6),
#     RHO_BOOT_Cpp_Parallel(x,399,7),
#     times = 100
# )
# RHO_BOOT_Cpp_Parallel(x,99,4)

*/

