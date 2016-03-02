
##' Adaptive t-test as described in Timmesfeld et al. (2007)
##'
##' Warning the current implementation seems to be numerically instable
##' @title Adaptive t-test
##' @param mpfr Whether to use high precision numbers from \code{Rmpfr}
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export
adaptive_ttest_os <- function(x,n1,n,ne,alpha=0.025,mpfr=FALSE) {
    if(n == ne){
        return(0.025 >= t.test(x,alternative='greater')$p.value)
    }
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    V1 <- sum(xs[[1]])
    U <- sum(xs[[1]]^2)
    tU <- sum(xs[[2]]^2)
    A <- clev(tU,U,V1,ne-n1,n,n1,alpha=alpha,mpfr=mpfr)
    A >= t.test(xs[[2]],alternative='greater')$p.value
}

##' Adaptive combination test of stage-wise t-tests using the inverse normal combination function. 
##'
##' @title Inverse normal adaptive t-test
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export
adaptive_invnorm_ttest_os <- function(x,n1,n,ne,alpha=0.025){
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    p1 <- t.test(xs[[1]],alternative='greater')$p.value
    p2 <- t.test(xs[[2]],alternative='greater')$p.value
    pnorm(sum(sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)),lower=FALSE) <= alpha 
}

##' Adaptive combination test of stage-wise t-tests using the inverse normal combination function. 
##'
##' @title Inverse normal adaptive wilcoxon test
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export
adaptive_invnorm_wilcoxtest_os <- function(x,n1,n,ne,alpha=0.025){
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    p1 <- wilcox.test(xs[[1]],alternative='greater')$p.value
    p2 <- wilcox.test(xs[[2]],alternative='greater')$p.value
    pnorm(sum(sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)),lower=FALSE) <= alpha 
}


##' Non-parametric combination of stage-wise test statistics. Combines stage-wise permutation p-values using some combination function; performs the test using the joint conditional permutation distribution of stage wise permutation p-values.
##'
##' @title adptive NPC test
##' @template onesample_sims
##' @param test_statistic Function that computes the test test statistic
##' @param combination_function Function to combine stage-wise (permutation) p-values
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
##' @export
adaptive_npcombtest_os <- function(x,n1,n,ne,test_statistic,combination_function=inverse_normal,perms=50000,alpha=0.025) {
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    gs <- split(sign(x)>0,rep(1:2,c(n1,ne-n1)))
    G <- omega(gs[[1]],gs[[2]],restricted=FALSE,B=perms)
    rB <- ncol(G)
    p1 <- 1-(rank(test_statistic(xs[[1]],G[1:n1,]))/(rB+1))
    p2 <- 1-(rank(test_statistic(xs[[2]],G[(n1+1):ne,]))/(rB+1))
    ct <- combination_function(p1,p2,n1/n,(n-n1)/n)
    sum(ct[1]>=ct)/length(ct) <= alpha
}










######## other adaptive tests
##################
##' Adaptive combination test of stage-wise 2 samples t-tests using the inverse normal combination function. 
##'
##' @title Inverse normal adaptive wilcoxon test
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export
adaptive_invnormtest_2s <- function(x,y,n1,n,ne,m1=n1,m=n,me=ne,alpha=0.025){
  xs <- split(x,rep(1:2,c(n1,ne-n1)))
  ys <- split(y,rep(1:2,c(m1,ne-m1)))
  p1 <- t.test(xs[[1]],ys[[1]],alternative='less')$p.value
  p2 <- t.test(xs[[2]],ys[[2]],alternative='less')$p.value
  alpha >= {sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)} %>% sum() %>% pnorm(lower=FALSE) 
}

##' Adaptive combination test of stage-wise 2 samples t-tests using the inverse normal combination function. 
##'
##' @title Inverse normal adaptive wilcoxon test
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export

adaptive_invnormtest_negbin_2s <- function(x,y,n1,n,ne,m1=n1,m=n,me=ne,alpha=0.025){
  xs <- split(x,rep(1:2,c(n1,ne-n1)))
  ys <- split(y,rep(1:2,c(m1,ne-m1)))
  sg1 <- summary(glm.nb(c(xs[[1]],ys[[1]])~rep(0:1,c(n1,m1))))
  sg2 <- summary(glm.nb(c(xs[[2]],ys[[2]])~rep(0:1,c(ne-n1,ne-m1))))
  p1 <- sg1$coefficients[2,"Pr(>|z|)"]
  s1 <- sg1$coefficients[2,"z value"]>0
  p2 <- sg2$coefficients[2,"Pr(>|z|)"]
  s2 <- sg2$coefficients[2,"z value"]>0
  p1 <- ifelse(s1,p1/2,1-p1/2)
  p2 <- ifelse(s2,p2/2,1-p2/2)
  alpha >= {sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)} %>% sum() %>% pnorm(lower=FALSE) 
}

##' Adaptive combination test of stage-wise 2 samples t-tests using the inverse normal combination function. 
##'
##' @title Inverse normal adaptive wilcoxon test
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export
adaptive_waldtest_2s <- function(x,y,n1,n,ne,m1=n1,m=n,me=ne,alpha=0.025){
  log(mean(y)/mean(x))/sqrt(1/sum(x)+1/sum(y))>qnorm(alpha,lower=F)
}
############################## other 2 samples

##' Adaptive combination test of stage-wise 2 samples t-tests using the inverse normal combination function. 
##'
##' @title Inverse normal adaptive wilcoxon test
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export

adaptive_invnorm_wilcoxtest_2s <- function(x,y,n1,n,ne,m1=n1,m=n,me=ne,alpha=0.025){
  xs <- split(x,rep(1:2,c(n1,ne-n1)))
  ys <- split(y,rep(1:2,c(m1,ne-m1)))
  p1 <- wilcox.test(xs[[1]],ys[[1]],alternative='less')$p.value
  p2 <- wilcox.test(xs[[2]],ys[[2]],alternative='less')$p.value
  alpha >= {sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)} %>% sum() %>% pnorm(lower=FALSE) 
}


