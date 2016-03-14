// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// all_combinations_cpp
IntegerMatrix all_combinations_cpp(int n, int k);
RcppExport SEXP adaperm_all_combinations_cpp(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    __result = Rcpp::wrap(all_combinations_cpp(n, k));
    return __result;
END_RCPP
}
// all_reassignments_cpp
LogicalMatrix all_reassignments_cpp(int n, int k);
RcppExport SEXP adaperm_all_reassignments_cpp(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    __result = Rcpp::wrap(all_reassignments_cpp(n, k));
    return __result;
END_RCPP
}
// diffmeanC
NumericVector diffmeanC(NumericVector x, IntegerMatrix g);
RcppExport SEXP adaperm_diffmeanC(SEXP xSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type g(gSEXP);
    __result = Rcpp::wrap(diffmeanC(x, g));
    return __result;
END_RCPP
}
// sumdiffC
NumericVector sumdiffC(NumericVector x, IntegerMatrix g);
RcppExport SEXP adaperm_sumdiffC(SEXP xSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type g(gSEXP);
    __result = Rcpp::wrap(sumdiffC(x, g));
    return __result;
END_RCPP
}
// meandiffC
NumericVector meandiffC(NumericVector x, IntegerMatrix g);
RcppExport SEXP adaperm_meandiffC(SEXP xSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type g(gSEXP);
    __result = Rcpp::wrap(meandiffC(x, g));
    return __result;
END_RCPP
}
// random_reassignments_cpp
IntegerMatrix random_reassignments_cpp(const IntegerVector g, const int nperm);
RcppExport SEXP adaperm_random_reassignments_cpp(SEXP gSEXP, SEXP npermSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const IntegerVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< const int >::type nperm(npermSEXP);
    __result = Rcpp::wrap(random_reassignments_cpp(g, nperm));
    return __result;
END_RCPP
}
// random_samples_cpp
NumericMatrix random_samples_cpp(const NumericVector xs, const int k, const int nsam);
RcppExport SEXP adaperm_random_samples_cpp(SEXP xsSEXP, SEXP kSEXP, SEXP nsamSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type nsam(nsamSEXP);
    __result = Rcpp::wrap(random_samples_cpp(xs, k, nsam));
    return __result;
END_RCPP
}
