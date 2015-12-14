// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// diffmeanC
NumericVector diffmeanC(NumericVector x, IntegerMatrix g);
RcppExport SEXP resamplingMCP_diffmeanC(SEXP xSEXP, SEXP gSEXP) {
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
RcppExport SEXP resamplingMCP_sumdiffC(SEXP xSEXP, SEXP gSEXP) {
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
RcppExport SEXP resamplingMCP_meandiffC(SEXP xSEXP, SEXP gSEXP) {
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
RcppExport SEXP resamplingMCP_random_reassignments_cpp(SEXP gSEXP, SEXP npermSEXP) {
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
RcppExport SEXP resamplingMCP_random_samples_cpp(SEXP xsSEXP, SEXP kSEXP, SEXP nsamSEXP) {
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
// combinations_cpp
IntegerMatrix combinations_cpp(const int n, const int k);
RcppExport SEXP resamplingMCP_combinations_cpp(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    __result = Rcpp::wrap(combinations_cpp(n, k));
    return __result;
END_RCPP
}
// subsamples_cpp
NumericMatrix subsamples_cpp(NumericVector xs, const int k);
RcppExport SEXP resamplingMCP_subsamples_cpp(SEXP xsSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    __result = Rcpp::wrap(subsamples_cpp(xs, k));
    return __result;
END_RCPP
}
// bincombinations_cpp
IntegerMatrix bincombinations_cpp(const int p);
RcppExport SEXP resamplingMCP_bincombinations_cpp(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    __result = Rcpp::wrap(bincombinations_cpp(p));
    return __result;
END_RCPP
}
