#include <Rcpp.h>
using namespace Rcpp;


// LEGGI QUESTA REFERENCE PER PROGRAMMARE
// http://adv-r.had.co.nz/Rcpp.html#rcpp-sugar
//
// ES. le funzioni di Sugar permettono di usare
//     operazioni vettoriali


// [[Rcpp::export]]
NumericVector Xresid(NumericVector X, NumericVector alphas) {
    unsigned int n = X.length();
    unsigned int order = alphas.length();
    NumericVector resid(n-order);

    for(unsigned int t = order; t < n; t++){
        NumericVector Xtmp = X[Rcpp::Range(t-order,t-1)];
        std::reverse(Xtmp.begin(), Xtmp.end());

        // std::cout << "t" << std::endl;
        // std::cout << t << std::endl;
        // std::cout << std::inner_product(Xtmp.begin(), Xtmp.end(), alphas.begin(), 0.0) << std::endl;
        resid[t - order] = X[t] - std::inner_product(Xtmp.begin(), Xtmp.end(), alphas.begin(), 0.0);
    }
    return resid;
}


// [[Rcpp::export]]
Rcpp::List Xmoments(NumericVector X, NumericVector alphas) {
    // unsigned int n = X.length();
    // unsigned int order = alphas.length();

    double mX = mean(X);
    double vX = var(X);

    NumericVector resid = Xresid(X, alphas);

    double mINN = (1 - sum(alphas))*mX; // (1 - sum(a))*mean(X)
    double vINN = var(resid);

    return List::create(
        _["meanX"]  = mX,
        _["varX"]  = vX,
        _["meanINN"]  = mINN,
        _["varINN"]  = vINN,
        _["resid"]  = resid
    );
}

/*** R
# aa <- c(0.1,0.2)
# asd <- Xmoments(rpois(500,2),aa)
# asd
# # mean(asd$resid);asd$meanINN
# plot(asd$resid,type = "l")
# abline(h=0,col="red")
#
# es INAR(p)
# x <- 1:10
# p <- 3
# aa <- c(0.2,0.4,0.1)
# Xresid(x,aa)
# for(t in (p+1):10){
#     i <- x[t] - x[(t-1):(t-p)]%*%aa
#     cat(i,"\n")
# }

*/
