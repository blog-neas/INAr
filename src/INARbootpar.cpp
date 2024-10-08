
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
