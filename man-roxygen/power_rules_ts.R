##' @param x1 First stage observations
##' @param y1 treatment group obseravtions
##' @param m Target mean value
##' @param target Target power
##' @param alpha Significance level
##' @param maxN Maximum sample size
##' @param rob_var Use robust variance estimate
##' @return Total sample size required to achieve the target power 
##' @details In case the robust variance estimate is used we compute the variance as \code{(iqr(c(x1-median(x1),y1-median(y1)))/1.349)^2} this leads to a slight underestimation of the required sample size - for the moment you could increase the target power to correct
