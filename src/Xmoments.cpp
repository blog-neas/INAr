#include <Rcpp.h>
using namespace Rcpp;


// LEGGI QUESTA REFERENCE PER PROGRAMMARE
// http://adv-r.had.co.nz/Rcpp.html#rcpp-sugar
//
// ES. le funzioni di Sugar permettono di usare
//     operazioni vettoriali


// // [[Rcpp::export]]
// NumericVector Xresid(NumericVector X, NumericVector alphas) {
//     unsigned int n = X.length();
//     unsigned int order = alphas.length();
//     NumericVector resid(n-order);
//
//     for(unsigned int t = order; t < n; t++){
//         NumericVector Xtmp = X[Rcpp::Range(t-order,t-1)];
//         std::reverse(Xtmp.begin(), Xtmp.end());
//
//         // std::cout << "t" << std::endl;
//         // std::cout << t << std::endl;
//         // std::cout << std::inner_product(Xtmp.begin(), Xtmp.end(), alphas.begin(), 0.0) << std::endl;
//         resid[t - order] = X[t] - std::inner_product(Xtmp.begin(), Xtmp.end(), alphas.begin(), 0.0);
//     }
//     return resid;
// }


// [[Rcpp::export]]
Rcpp::List Xresid(NumericVector X, NumericVector alphas, double mINN, double vINN) {
    unsigned int n = X.length();
    unsigned int order = alphas.length();

    NumericVector alphasmix = alphas*(1-alphas);
    NumericVector resid(n-order);
    NumericVector stdresid(n-order);

    double innProdAplhas;
    double innProdAplhasmix;
    double mXcond;
    double vXcond;

    for(unsigned int t = order; t < n; t++){
        NumericVector Xtmp = X[Rcpp::Range(t-order,t-1)];
        std::reverse(Xtmp.begin(), Xtmp.end());

        // std::cout << "t" << std::endl;
        // std::cout << t << std::endl;
        // std::cout << std::inner_product(Xtmp.begin(), Xtmp.end(), alphas.begin(), 0.0) << std::endl;
        innProdAplhas = std::inner_product(Xtmp.begin(), Xtmp.end(), alphas.begin(), 0.0);
        innProdAplhasmix = std::inner_product(Xtmp.begin(), Xtmp.end(), alphasmix.begin(), 0.0);
        mXcond = innProdAplhas + mINN;
        vXcond = innProdAplhasmix + vINN;

        resid[t - order] = X[t] - innProdAplhas;
        stdresid[t - order] = (X[t] - mXcond)/sqrt(vXcond);

        // std::cout << "t" << std::endl;
        // std::cout << resid[t - order] << std::endl;
        // std::cout << "std" << std::endl;
        // std::cout << stdresid[t - order] << std::endl;
    }

    return  List::create(
        _["res"] = resid,
        _["stdres"] = stdresid
    );
}


// [[Rcpp::export]]
Rcpp::List Xmoments(NumericVector X, NumericVector alphas) {
    // unsigned int n = X.length();
    // unsigned int order = alphas.length();

    double mX = mean(X);
    double vX = var(X);

    double mINN = (1 - sum(alphas))*mX; // (1 - sum(a))*mean(X)
    double vINN = vX*(1 - sum(pow(alphas,2))) - mX*sum(alphas*(1-alphas));

    List RES = Xresid(X, alphas, mINN, vINN);
    NumericVector resid = RES["res"];
    NumericVector stdresid = RES["stdres"];

    // questi indici mi sa che non hanno molto senso.... o no?
    // double mINN2 = mean(resid);
    // double vINN2 = var(resid);

    return List::create(
        _["meanX"]  = mX,
        _["varX"]  = vX,
        _["meanINN"]  = mINN,
        _["varINN"]  = vINN,
        _["resid"]  = resid,
        _["stdresid"]  = stdresid
    );
}

/*** R
# aa <- c(0.1,0.2,0.4)
# asd <- Xmoments(rpois(10000,2),aa)
# asd$meanX;asd$varX
# asd$meanINN;asd$varINN;
# # asd$meanINN2;asd$varINN2
# plot(asd$resid,type = "l")
# abline(h=0,col="red")
# plot(asd$stdresid,type = "l")
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
