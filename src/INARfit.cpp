#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
