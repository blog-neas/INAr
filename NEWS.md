# INAr 0.2.3

* Added  the `disc_unif`, `mix_bin` and `mix_bin_negbin` arrival options into the `genINAR`function.

* Applying some modifications to the `SMC_Cpp`function for the Generalized Poisson case. In particular, I am trying to figure out how to deal with some limit cases in which the test statistic is not defined, that is when $p_{x{t}-1}$ and/or $p_{x{t}}$ are zero.

* [testing] Playing with some solutions for parallel computing in the C++ code.


# INAr 0.2.2

* Added Rho test statistic that will be used as benchmark for the other tests.

* Added a `INARtest` wrapper function to call the INAR tests.


# INAr 0.2.1

* [testing] Added Harris-McCabe test statistic.
    * Bootstrap procedures are still on development phase.


# INAr 0.2.0

* Added parametric and semiparametric bootstrap procedures for the Sun-McCabe Score test.


# INAr 0.1.1

* [testing] Preparing the support of the Harris-McCabe test statistics for the upcoming 0.2 update


# INAr 0.1.0

* Added parametric and semiparametric bootstrap procedures for the Sun-McCabe Score test in case of Generalized Poisson innovations.


# INAr 0.0.12

* Substituted the `NMF` (function `fcnnls`) library with `RcppML` (function `nnls`), now the package Biobase is no more an upstream dependency ([issue 1](https://github.com/blog-neas/INAr/issues/1)).

* Solved a bug in `INARp_cpp`: alphas and lagged values were inverted.

* Added two new blocks of scripts: `genericfuns.R` and `utils.R` for the upcoming 0.1 update.

* Preparing the main frontend function `INARfit` for the upcoming 0.1 update:
  * Introduced the new `INAR` class;
  * Added some generic functions;
  * YW and CLS estimation are available;
  * Minor changes and few corrections.

* Added the Negative Binomial SMC parametric bootstrap test and the PIT experimental.
  * NOTE: even if they seem to work properly, SMC tests need some cleaning and a thorough check!


# INAr 0.0.11

* Added a new test in `test.R` that follow the same concept of the previous one:
    * Added the `SMCboot.test` function that computes the semiparametric or parametric bootstrap Sun-McCabe Score test statistics (with Poisson or Negative Binomial arrivals for the moment). The function returns an object of class `htest`.
    * In future the the C++ routines that compute the tests will become internal.


# INAr 0.0.10

* Added a new script, namely `test.R`, that will include all the front-end test functions.
    * Added the `SMC.test` function that computes the Sun-McCabe Score test statistics (with Poisson or Negative Binomial arrivals for the moment). The function returns an object of class `htest`.
    * In future the bootstrapped version of the above tests will be added and the C++ routines that compute the tests will become internal.


# INAr 0.0.9

* Preparing for the first CRAN submission:
    * Minor changes and few modifications;
    * Script cleaning.


# INAr 0.0.8

* Cleaning and few modifications.

* Changed the formula regarding he Negative Binomial parameters' estimation acording with Sun, McCabe (2013).

* Parametric bootstrap of INAR with Negative Binomial innovations is now available, although it needs some testing.


# INAr 0.0.7

* Added downloads dataset, source Weiss (2008).

* First build of the package vignette. 


# INAr 0.0.6

* Added package sticker.

* Implementing the Sun-McCabe bootstrap test. This cose is still in development and works properly only for INAR(1) processes. 


# INAr 0.0.5

* Improvement of the `INARfit.R` code to fit INAR(p) models. Now `INARfit()` performs a full Y-W estimation from a Poisson INAR(r) family, following the results of Du and Li. Some additional Rcpp utility functions (script `Xmoments.cpp`) have been added:
    * `Xmoments()` [in development], compute the first two moments for the original series and the residual series. As output it returns mean and variance of both the starting and residual series, and the estimated residual series;
    * `Xresid()` [in development], generates the series of residual values.


# INAr 0.0.4

* First `INARfit.R` code to fit INAR(p) models:
    * `INARfit()` [in development], fitting an INAR(p) process, by using several procedures. At the momemt is hardcoded and works only for the Poisson case and only YW is provided;
    * `est_mom()` [in development], estimation of innovations' parameters. At the momemt is hardcoded and works only for the Poisson case.


# INAr 0.0.3

* Generalization of genINAR function, now it generates INAR(p) models nstead of INAR(1).

* Added stationarity condition check in `genINAR()` function.

* The old `par` input vector contained both the thinning operator (at the first position) and innovations' parameters, now this vector is split is two: `a` and `par`, where:
    * `a` contains the p thinning parameters of the INAR(p) to be generated
    * `par` contains exclusively the innovations' parmeters

* Development of the C++ part to generate INAR(p) processes
    * deleted INAR1\_gen.cpp and included the routine `INAR1\_ cpp` in  INARp\_gen.cpp
    * development of the more general routine `INARp\_ cpp`; 
    * the line `sim = clone(resid)` was added to avoid the shallow copy effect;
    * INARp_gen.cpp contains the routine to generate an INAR(p) process;
    * at the moment `INAR1\_ cpp` si obsolete, after some testing it will be deleted.
 

# INAr 0.0.2

* Updated `README.md` file.

* Added a `NEWS.md` file to track main changes among different versions.

* Added references.

* `DESCRIPTION` file updated
    * additions: `URL`, `Roxygen`, `Depends: R (>= 4.2.1)`, `RdMacros`, `LinkingTo`, `RoxygenNote`;
    * modifications: `Depends: R (>= 4.2.1)`, `Imports: Rcpp (>= 1.0.0)`, `RcppArmadillo`, `MASS`, `Rdpack`.


# INAr 0.0.1

* Initialization of the package, first settings.

* First working version of the package:
    * use of `Rcpp` and `RccpARmadillo` to add C++ code;
    * a preliminary version of the function `inarGEN()` is implemented.
