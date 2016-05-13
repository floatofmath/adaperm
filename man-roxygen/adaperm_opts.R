##' @param x Observations
##' @param g Treatment assignments (if \code{NULL}) a one-sample test will be performed
##' @param n1 First stage control-group sample size
##' @param n Pre-planned overall sample size
##' @param m1 First stage treatment-group sample size (will be ignored if \code{g} is \code{NULL})
##' @param m Pre-planned overall treatment-group sample size (will be ignored if \code{g} is \code{NULL})
##' @param test_statistic Test statistic
##' @param alpha Significance level
##' @param cer_type what type of conditional error rate function should be used (see \code{\link{permutation_cer}}) for detailcer_type 
##' @param atest_type if 'CER' compute only conditional error rate, else type of adaptive test should be performed (see \code{\link{perm_test}} for details) 
##' @param permutations Number of permutations to use
##' @param stratified should permutation be stratified by stage
