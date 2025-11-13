#include <Rcpp.h>
using namespace Rcpp;

// estimate a series to compute variance of ZIpv test
// [[Rcpp::export]]
double series_varpv(double mu, double alpha, int max_j = 2000) {
    double sum = 0.0;
    double expon = 1.0; // to compute mu^j / j! incrementally
    // check if alpha <1
    if (alpha < 0 || alpha >= 1) {
        Rcout << "alpha must be in [0, 1)" << std::endl;
        return NA_REAL;
    }

    for (int j = 1; j <= max_j; ++j) {
        // Rcout << "j: " << j << std::endl;

        expon *= mu / j;  // updates mu^j / j!
        double numerator = expon * pow(alpha, j);
        double denominator = 1.0 - pow(alpha, j);
        sum += numerator / denominator;

        // optional: stop early if the sum becomes negligible
        if (fabs(numerator / denominator) < 1e-15) break;
    }

    return sum;
}

/*** R
# series_varpv(2, 0.5, 1000)
# series_varpv(2, 0, 1000)
# series_varpv(2, 0.99, 1000)
# series_varpv(2, 1, 1000)
*/
