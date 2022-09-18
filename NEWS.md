
# INAr 0.0.3

* Generalization of genINAR function, now it generates INAR(p) models nstead of INAR(1).
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
