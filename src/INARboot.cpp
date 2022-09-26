// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

using namespace R;
using namespace Rcpp;

//
// NB. HO INCLUSO ANCHE LO SCRIPT PRESENTE IN statistica_SunMcCabe.cpp
//

//
// !!!! DA IMPLEMENTARE ANCORA !!!!
//

// [[Rcpp::export]]
NumericVector sunMC_Cpp(NumericVector x, unsigned int method){
  int n = x.length();

  NumericVector out(2);
  Range idx = seq(1,n-1);
  Range idx_1 = seq(0,n-2);

  NumericVector  xsum = x[idx];
  NumericVector  xsum_1 = x[idx_1];

  // double var_x = var(noNA(x));

  // lambda_hat <- mean(x)
  //
  //   mu_g <- NA # nun ce serv
  //   sd_g <- NA # nun ce serv
  //
  //   xsumL <- xsum-lambda_hat;
  // xsumL_1 <- xsum_1-lambda_hat;
  // NUM <- xsumL%*%xsumL_1
  //
  //   stat <- NUM/lambda_hat
  //   stat <- (n^(-1/2))*stat
  //
  //   pval <- 1-pnorm(stat)

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
    // printf("1) %f, %f %f \n",mu_x,var_x,std::abs(var_x-mu_x));
    // printf("2) %f, %f \n",pow(mu_x,2),var_x-mu_x);

    double r_hat = pow(mu_x,2)/std::fabs(var_x-mu_x); // # trick
    double p_hat = std::fabs(var_x - mu_x)/var_x; // # trick
    // printf("%f, %f \n",r_hat,p_hat);

    // NumericVector ind_val (n-1); // n-1 = xsum.length()
    // LogicalVector NoZero = xsum > 0;
    // NumericVector xsumNoZ = xsum[NoZero];
    // unsigned int n_v = xsumNoZ.length();
    // printf("%i \n",n_v);
    // NumericVector r_vec(n_v,r_hat);
    // NumericVector p_vec(n_v,p_hat);
    // printf("%f \n",r_vec[4]);
    // printf("%f \n",p_vec[7]);
    //
    // printf("%f, %f \n",r_hat,p_hat);
    // ind_val[NoZero] = xsumNoZ/((r_vec+xsumNoZ-1)*p_hat);
    // NumericVector asd = xsumNoZ/((r_hat+xsumNoZ-1)*p_hat);

    NumericVector ind_val = ifelse(xsum > 0, xsum/((r_hat+xsum-1)*p_hat), 0);

    // for(int i=0; i<asd.length(); ++i){
    //   Rprintf("the value of v[%i] : %f | xsum[%i] : %f  | (%f) \n", i, asd[i], i, xsum[i],xsum[i]/((r_hat+xsum[i]-1)*p_hat));
    // }

    //
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



// [[Rcpp::export]]
NumericVector sunMC_semiparBOOT_Cpp(NumericVector x, int B, unsigned int method){
  int n = x.length();
  // NumericVector out(3);
  // NumericVector MB(B);
  // NumericVector VB(B);

  NumericVector s_temp(B);
  for(int i = 0; i < B; i++){

    NumericVector xb(n);
    // print(xb);

    // check!
    LogicalVector id(n);
    id = xb==xb[0];

    while(Rcpp::all(id).is_true()) {
      xb = RcppArmadillo::sample(x,n,true);
      id = xb==xb[0];
    }

    s_temp[i] = sunMC_Cpp(xb,method)[0];
  }

  return s_temp;
}

//
// !!!! DA IMPLEMENTARE !!!!
//
// [[Rcpp::export]]
NumericVector sunMC_parBOOT_Cpp(NumericVector x, int B, unsigned int method){
    int n = x.length();
    double lambda_x = mean(noNA(x)); // = lambda

    // NumericVector out(3);
    // NumericVector MB(B);
    // NumericVector VB(B);

    NumericVector s_temp(B);
    for(int i = 0; i < B; i++){

        NumericVector xb(n);
        // print(xb);

        // check!
        LogicalVector id(n);
        id = xb==xb[0];

        while(Rcpp::all(id).is_true()) {
            // POISSON
            xb = Rcpp::rpois(n,lambda_x);
            id = xb==xb[0];
        }

        s_temp[i] = sunMC_Cpp(xb,method)[0];
    }

    return s_temp;
}


// [[Rcpp::export]]
List sunMCtest_boot(NumericVector X, unsigned int arrival, unsigned int method, int B){
    // int n = X.length();

    double Smc = sunMC_Cpp(X,arrival)[0];
    double B_stat;
    double B_pval;
    NumericVector SmcB(B);

    if(method == 1){
        // SEMIPARAMETRIC
        SmcB = sunMC_semiparBOOT_Cpp(X,B,arrival);
    }
    else if(method == 2){
        // PARAMETRIC
        SmcB = sunMC_parBOOT_Cpp(X,B,arrival);
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
# S <- sunMC_Cpp(x,method)
# Sb <- sunMC_semiparBOOT_Cpp(x,B,method)
# mean(S[1] > Sb) # pval boot
# S[2] # pval normale
# # test 2: FUNZIONE FINALE WRAPPED
# arrival <- 1
# sunMCtest_boot(x,arrival,method,B)
*/