# INAr 0.3.4

* [LP] Updating estimation procedures. In particular, I am working on the extension of the Yule-Walker(YW) and Conditional Least Squares (CLS) estimation procedure from only INAR(1) to more general INAR(p) processes. Now the two procedures are available for the Poisson, Negative Binomial and Generalized Poisson cases, they will be extended to the Katz case in the next updates.

* [LP] Some minor bug fixes in the estimation procedures and in the C++ code for the generation of INAR(p) processes.

* [LP] Writing the main vignette to show how to use the package and to illustrate the main features of the INAR processes. The vignette is still in development phase, but it will be available in the next updates.

* [LP] Adding ggacf() / ggpacf() ggplot-based ACF/PACF helpers.

* [LP] Updating package metadata (DESCRIPTION/NEWS/NAMESPACE) for the new release and vignette dependencies.

* [LP] Updating github actions.


# INAr 0.3.3

* [MTR] Added some vignettes *in development phase* to show how to use the package and to illustrate the main features of the INAR processes.

* [MTR] Added the cryptosporidiosis infection count data.


# INAr 0.3.2

* [LP] Fixed minor bugs in the test functions, some tests were not showing the alternative hypothesis in the output.

* [LP] Implemented some *experimental* functions to work with multidimensional INAR processes, in particular:
    * `genMINAR` (in MINAR.R): Frontend function to generate MINAR(p) models.
    * `estimMYW` (in MINAR.R): Function to fit MINAR(p) models by using the Yule-Walker procedure.
    * `MINARp_gen_cpp` (in MINARp_gen.cpp): C++ routine to generate MINAR(p) processes.

* [LP] Initial development of the Katz distribution computation, in particular:
    * `dKatz` (in Katz.R): Function to compute the Katz distribution probability mass function.
    * `pKatz` (in Katz.R): Function to compute the Katz distribution cumulative distribution function.
    * `qKatz` (in Katz.R): Function to compute the Katz distribution quantile function.
    * `rKatz` (in Katz.R): Function to generate random numbers from the Katz distribution.

# INAr 0.3.1

* [LP] Implemented the Harris-McCabe bootstrap test.
    * `HMCtest`: Main frontend function to perform the Harris-McCabe test.

* [LP] Fixed some bugs in the C++ code for the bootstrap procedures.

* [LP] Cleaning of the overall scripts from comments and debugging lines.

* [LP] Implemented a generic function `print.INARtest` to print the results of all the test functions.


# INAr 0.3.0

* [LP] Overhaul of the package structure and code organization. The package is now organized into multiple R scripts, each dedicated to specific functionalities such as data generation, model fitting, statistical tests and utility functions. List of the main functions:
    * `genINAR`: Main frontend function to generate INAR(p) models.
    * `INAR`: Main frontend function to fit INAR(p) models.
    * `DItest`: Main frontend function to perform the Dispersion Index test.
    * `ZItest`: Main frontend function to perform the Zero-Inflation test.
    * `ZIDItest`: Main frontend function to perform the combined Zero-Inflation and Dispersion Index test.
    * `SMCtest`: Main frontend function to perform the Sun-McCabe Score test.
    * `HMCtest`: Main frontend function to perform the Harris-McCabe test (under development).

* [LP] New function organization to improve scalability and maintainability. The function dependency structure is described in the readme file.

* [LP] Reshaping all the C++ code to improve efficiency and speed.

* [LP] Parallel C++ scripts for bootstrap procedures under development phase.


# INAr 0.2.3

* [LP] Removed the the `good` package since it is not available anymore on CRAN. In next updates I'm considering to substitute the `good::rgood` function by using the `VGAM` package.

* [LP] Added  the `disc_unif`, `mix_bin` and `mix_bin_negbin` arrival options into the `genINAR`function. Solving minor bugs.

* [LP] General debugging and modifications only for the Generalized Poisson case in `SMC_Cpp` function: In particular, I am trying to figure out how to deal with some limit cases in which the test statistic is not defined, that is when $p_{x{t}-1}$ and/or $p_{x{t}}$ are zero.

* [LP] Introduced the `SMCtest` wrapper function for the Sun-McCabe Score test. The `INARtest` function will be deprecated in the next few updates.

* [LP] Playing with some potential solutions for parallel computing in the C++ code.


# INAr 0.2.2

* [LP] Added Rho test statistic that will be used as benchmark for the other tests.

* [LP] Added a `INARtest` wrapper function to call the INAR tests.


# INAr 0.2.1

* [LP] Added Harris-McCabe test statistic (under development).
    * Bootstrap procedures are still on development phase.


# INAr 0.2.0

* [LP] Added parametric and semiparametric bootstrap procedures for the Sun-McCabe Score test.


# INAr 0.1.1

* [LP] Testing the support of the Harris-McCabe test statistics for the upcoming 0.2 update


# INAr 0.1.0

* [LP] Added parametric and semiparametric bootstrap procedures for the Sun-McCabe Score test in case of Generalized Poisson innovations.


# INAr 0.0.12

* [LP] Substituted the `NMF` (function `fcnnls`) library with `RcppML` (function `nnls`), now the package Biobase is no more an upstream dependency ([issue 1](https://github.com/blog-neas/INAr/issues/1)).

* [LP] Solved a bug in `INARp_cpp`: alphas and lagged values were inverted.

* [LP] Added two new blocks of scripts: `genericfuns.R` and `utils.R` for the upcoming 0.1 update.

* [LP] Preparing the main frontend function `INARfit` for the upcoming 0.1 update:
  * Introduced the new `INAR` class;
  * Added some generic functions;
  * YW and CLS estimation are available;
  * Minor changes and few corrections.

* [LP] Added the Negative Binomial SMC parametric bootstrap test and the PIT experimental.
  * NOTE: even if they seem to work properly, SMC tests need some cleaning and a thorough check!


# INAr 0.0.11

* [LP] Added a new test in `test.R` that follow the same concept of the previous one:
    * Added the `SMCboot.test` function that computes the semiparametric or parametric bootstrap Sun-McCabe Score test statistics (with Poisson or Negative Binomial arrivals for the moment). The function returns an object of class `htest`.
    * In future the the C++ routines that compute the tests will become internal.


# INAr 0.0.10

* [LP] Added a new script, namely `test.R`, that will include all the front-end test functions.
    * Added the `SMC.test` function that computes the Sun-McCabe Score test statistics (with Poisson or Negative Binomial arrivals for the moment). The function returns an object of class `htest`.
    * In future the bootstrapped version of the above tests will be added and the C++ routines that compute the tests will become internal.


# INAr 0.0.9

* [LP] Preparing for the first CRAN submission:
    * Minor changes and few modifications;
    * Script cleaning.


# INAr 0.0.8

* [LP] Cleaning and few modifications.

* [LP] Changed the formula regarding he Negative Binomial parameters' estimation acording with Sun, McCabe (2013).

* [LP] Parametric bootstrap of INAR with Negative Binomial innovations is now available, although it needs some testing.


# INAr 0.0.7

* [LP] Added downloads dataset, source Weiss (2008).

* [LP] First build of the package vignette. 


# INAr 0.0.6

* [LP] Added package sticker.

* [LP] Implementing the Sun-McCabe bootstrap test. This cose is still in development and works properly only for INAR(1) processes. 


# INAr 0.0.5

* [LP] Improvement of the `INARfit.R` code to fit INAR(p) models. Now `INARfit()` performs a full Y-W estimation from a Poisson INAR(r) family, following the results of Du and Li. Some additional Rcpp utility functions (script `Xmoments.cpp`) have been added:
    * `Xmoments()` [in development], compute the first two moments for the original series and the residual series. As output it returns mean and variance of both the starting and residual series, and the estimated residual series;
    * `Xresid()` [in development], generates the series of residual values.


# INAr 0.0.4

* [LP] First `INARfit.R` code to fit INAR(p) models:
    * `INARfit()` [in development], fitting an INAR(p) process, by using several procedures. At the momemt is hardcoded and works only for the Poisson case and only YW is provided;
    * `est_mom()` [in development], estimation of innovations' parameters. At the momemt is hardcoded and works only for the Poisson case.


# INAr 0.0.3

* [LP] Generalization of genINAR function, now it generates INAR(p) models nstead of INAR(1).

* [LP] Added stationarity condition check in `genINAR()` function.

* [LP] The old `par` input vector contained both the thinning operator (at the first position) and innovations' parameters, now this vector is split is two: `a` and `par`, where:
    * `a` contains the p thinning parameters of the INAR(p) to be generated
    * `par` contains exclusively the innovations' parmeters

* [LP] Development of the C++ part to generate INAR(p) processes
    * deleted INAR1\_gen.cpp and included the routine `INAR1\_ cpp` in  INARp\_gen.cpp
    * development of the more general routine `INARp\_ cpp`; 
    * the line `sim = clone(resid)` was added to avoid the shallow copy effect;
    * INARp_gen.cpp contains the routine to generate an INAR(p) process;
    * at the moment `INAR1\_ cpp` si obsolete, after some testing it will be deleted.
 

# INAr 0.0.2

* [LP] Updated `README.md` file.

* [LP] Added a `NEWS.md` file to track main changes among different versions.

* [LP] Added references.

* [LP] `DESCRIPTION` file updated
    * additions: `URL`, `Roxygen`, `Depends: R (>= 4.2.1)`, `RdMacros`, `LinkingTo`, `RoxygenNote`;
    * modifications: `Depends: R (>= 4.2.1)`, `Imports: Rcpp (>= 1.0.0)`, `RcppArmadillo`, `MASS`, `Rdpack`.


# INAr 0.0.1

* [LP] Initialization of the package, first settings.

* [LP] First working version of the package:
    * use of `Rcpp` and `RccpARmadillo` to add C++ code;
    * a preliminary version of the function `inarGEN()` is implemented.

### Legend 

* [LP]: Lucio Palazzo

* [MTR]: Mariateresa Russo

