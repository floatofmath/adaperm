##' Split sample of a single group experiment into stages
##'
##' @title Split single group sample
##' @param x Observations
##' @param n1 First stage sample size
##' @param n Pre-planned total sample size
##' @return list with stagewise observations and stagewise sign-indicators
##' @author Florian Klinglmueller
split_sample_os <- function(x,n1,n){
    g <- sign(x)
    x <- abs(x)
    ne <- length(x)
    if(ne>n){
        xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
        gs <- split(g>0,rep(1:3,c(n1,n-n1,ne-n)))
        return(list(xs,gs))
    } else {
        xs <- split(x,rep(1:2,c(n1,ne-n1)))
        xs$`3` <- numeric(0)
        gs <- split(g>0,rep(1:2,c(n1,ne-n1)))
        gs$`3` <- logical(0)
        return(list(xs,gs))
    }
}

##' Split sample of a two group experiment into stages
##'
##' @title Split two-group sample
##' @param x Observations control group
##' @param y Observations treatment group
##' @param n1 First stage sample size
##' @param n Pre-planned total sample size
##' @param m1 First stage sample size
##' @param m Pre-planned total sample size
##' @return list with stagewise observations and stagewise sign-indicators
##' @author Florian Klinglmueller
split_sample_ts <- function(x,y,n1,n,m1,m){
    ne <- length(x)
    me <- length(y)
    if(ne>n){
        xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
    } else {
        if(me>m) stop('Stages with controls only not supported')
        xs <- split(x,rep(1:2,c(n1,ne-n1)))
        xs[[3]] <- numeric(0)
    }
    if(me>m){
        ys <- split(y,rep(1:3,c(m1,m-m1,me-m)))
    } else {
        if(ne>n) stop('Stages with treatments only not supported')
        ys <- split(y,rep(1:2,c(m1,me-m1)))
        ys[[3]] <- numeric(0)
    }
    gs <- lapply(1:length(xs),function(i) rep(0:1,c(length(xs[[i]]),length(ys[[i]]))))
    xs <- lapply(1:length(xs),function(i) c(xs[[i]],ys[[i]]))
    return(list(xs,gs))
}

##' Make contrast matrix for all (pairwise) many-to-one comparisons
##'
##' @param g Vector of treatment assignments
##' @title Contrasts for many-to-one comparisons
##' @export
##' @author Florian Klinglmueller
make_pw_contrasts <- function(g,control=0){
    levels <- unique(g)
    sapply(levels[levels!=control],function(l) (g == l) - (g == control))
}

##'  Simulate normal heteroscedastic data with scaling of the variance drawn from a chi-square distribution with \
##'
##' @title Generate heteroskedastic data
##' @param n  number of observations 
##' @param df degrees of freedom for the variance distribution
##' @param delta non-centrality parameter
##' @export
rnorm_hetero <- function(n,df=1,ncp=0){
  rnorm(n)/sqrt(rchisq(n,df)/df) + ncp
}


##' Simulate from a contaminated normal distribution. A proportion \code{csd} of random samples are scaled by \code{csd} and randomly shifted to the left or right by \code{cshift}.  
##'
##' @title Generate contaminated normally distributed data
##' @param n number of observations
##' @param mean mean value
##' @param sd standard deviation 
##' @param cprop proportion of contaminated samples
##' @param cshift mean shift of contaminated samples
##' @param csd standard deviation of contaminated samples
##' @return vector with \code{n} pseudo random numbers
##' @author Florian Klinglmueller
##' @export
rnorm_cont <- function(n,mean=0,sd=1,cprop=.1,cshift=3,csd=sd){
    cont <- sample(-1:1,n,prob=c(cprop/2,1-cprop,cprop/2),replace=TRUE)
    out <- rnorm(n,sd=sd) * (csd/sd)^abs(cont) + cshift*cont + mean
    out
}

##' Template function to compute the Type 1 error of a trial design
##'
##' @title Type 1 error 
##' @param X Random number generater for first stage observations
##' @param X2 Random number generator for second stage observations (if NULL first stage generator will be used)
##' @param n pre-planned total sample size
##' @param n1 first stagen sample size
##' @param B number of permutations to be used to compute the permutation null distribution
##' @param MCMC number of simulation runs for Type I error estimation
##' @param one_sample should a one-sample test be simulated
##' @param ... additional functions to conditional error rate function
##' @return matrix of p-values from different test procedures
##' @author Florian Klinglmueller
##'
##' @export
type1 <- function(X=rnorm,X2=NULL,n,n1,stat=meandiff,B,MCMC,one_sample=FALSE,...){
    x <- X(n)
    restricted=TRUE
    if(one_sample){
        x <- abs(x)
    }
    pm <- dplyr::bind_rows(lapply(1:MCMC,function(i) {
                                      g <- sample(rep(c(-1,1),each=n1/2))
                                      y <- if(is.null(X2)) x else c(x[1:n1],X2(n-n1))
                                      if(one_sample) y <- abs(y)
                                      list(preplanned=permutation_CER(x[1:n1],g,x[(n1+1):n],B=B,stat=stat,restricted=!one_sample),
                                           influenced=permutation_CER(y[1:n1],g,y[(n1+1):n],B=B,stat=stat,restricted=!one_sample),
                                           normal=normal_CER(x[1:n1],g,n,one_sample=one_sample,...),
                                           ttest=t_CER(x[1:n1],g,n,one_sample=one_sample,...))
        }))
    pm
}


##' Template function to compute the Type 1 error of a trial design
##'
##' @title Type 1 error 
##' @param X Random number generater for first stage observations
##' @param X2 Random number generator for second stage observations (if NULL first stage generator will be used)
##' @param ncp Non-centrality (shift) paramater 
##' @param n pre-planned total sample size
##' @param n1 first stagen sample size
##' @param stat 
##' @param B number of permutations to be used to compute the permutation null distribution
##' @param MCMC number of simulation runs for Type I error estimation
##' @param one_sample should a one-sample test be simulated
##' @param ... additional functions to conditional error rate function
##' @return matrix of p-values from different test procedures
##' @author Florian Klinglmueller
##'
##' @export
compare_power <- function(X=rnorm,X2=NULL,ncp,n,n1,stat=meandiff,B,MCMC,one_sample=FALSE,...){
    x <- X(n)
    if(one_sample) x <- abs(x)
    pm <- dplyr::bind_rows(lapply(1:MCMC,function(i) {
                                      g <- sample(rep(c(-1,1),each=n1/2))
                                      y <- if(is.null(X2)) x else c(x[1:n1],X2(n-n1))
                                      if(one_sample) y <- abs(y)
                                      list(preplanned=permutation_CER(x[1:n1],g,x[(n1+1):n],B=B,stat=stat,restricted=!one_sample),
                                           influenced=permutation_CER(y[1:n1],g,y[(n1+1):n],B=B,stat=stat,restricted=!one_sample),
                                           normal=normal_CER(x[1:n1],g,n,one_sample=one_sample,...),
                                           ttest=t_CER(x[1:n1],g,n,one_sample=one_sample,...))
                                  }))
    pm
}

##' Simulate trial data
##'
##' @title Simulate trials
##' @param MCMC Number of trials to simulate
##' @param n1 First stage sample size
##' @param n Total preplanned sample size
##' @param rfs First stage random number generator
##' @param rss Second stage random number generator
##' @param ncp Non-centrality parameter
##' @param n3 Third stage sample size (maybe a function of \code{x},\code{g2})
##' @param cond Condition on first and second stage observations ("both"), only on first stage observations ("first"), or unconditional of observations
##' @return list object with simulated data
##' @export
##' @author Florian Klinglmueller
simulate_trials <- function(MCMC,n1,n,n3=0,r=1/2,rfs=rnorm,rss=rnorm,ncp=0,cond=c("both","first","none"),restricted=TRUE){
    cond <- match.arg(cond)[1]
    n2 <- n-n1
    ns <- c(n1,n2)
    if(!isTRUE(all.equal(ns*r,as.integer(ns*r),check.attributes=FALSE))) { warning(paste("Relative group size",r,"do not permit whole numbered group sizes:",paste(ns*r,collapse=' ')," will be rounded!")) }
    ks <- floor(ns*r)
    if(ncp != 0 && cond != "none") { stop("Can not condition on observations under the alternative") }
    if(cond == "both"){
        x1 <- rfs(n1)
        x2 <- rss(n2)
        ans <- list(x=c(x1,x2),g=strat_reassignments(ns,ks,restricted=restricted,B=MCMC))
    } else if(cond == "first") {
        g1 <- strat_reassignments(ns[1],ks[1],restricted=restricted,B=MCMC)
        g2 <- if(restricted) matrix(c(rep(0L,ns[2]-ks[2]),rep(1L,ks[2])),ncol=ncol(g1),nrow=ns[2]) else matrix(sample(0L:1L,ns[2]*ncol(g1),replace=TRUE,prob=c(1-r,r)),nrow=ns[2])
        x1 <- rfs(n1)
        x2 <- matrix(rss(n2*ncol(g1)),nrow=n2)
        ans <- list(x=rbind(x1,x2),
                    g=rbind(g1,g2))
    } else if(cond == "none") {
        g1 <- if(restricted) matrix(c(rep(0L,ns[1]-ks[1]),rep(1L,ks[2])),nrow=ns[1],ncol=MCMC) else matrix(sample(0L:1L,ns[1]*MCMC,replace=TRUE,prob=c(1-r,r)),nrow=ns[1])
        g2 <- if(restricted) matrix(c(rep(0L,ns[2]-ks[2]),rep(1L,ks[2])),nrow=ns[2],ncol=MCMC) else matrix(sample(0L:1L,ns[2]*MCMC,replace=TRUE,prob=c(1-r,r)),nrow=ns[2])
        x1 <- matrix(rfs(n1*MCMC),nrow=n1)
        x2 <- matrix(rss(n2*MCMC),nrow=n2)
        x1 <- x1 + ncp*g1
        x2 <- x2 + ncp*g2
        ans <- list(x=rbind(x1,x2),
                    g=rbind(g1,g2))
        
    }
    if(is.function(n3)){
        n3 <- n3(x1,ans$g[,1:n1])
        k3 <- floor(n3*r)
        ans  <- within(ans,
                       {
                           if(restricted){
                               g3 <- lapply(1:ncol(g),function(i) c(rep(0L,n3[i]-k3[i])),rep(1L,k3))
                           } else {
                               g3 <- lapply(1:ncol(g),function(i) sample(0L:1L,n3[i],replace=TRUE,prob=c(1-r,r)))
                           }
                           x3 <- lapply(1:ncol(g),rss(n3[i])+ncp*g[i])
                       })
    } else if(n3>0) {
        k3 <- floor(n3*r)
        ans <- within(ans,
                      {
                          x3 <- matrix(rss(n3*ncol(g)),nrow=n3)+ncp*g3
                          g3 <- if(restricted) matrix(c(rep(0L,n3-k3),rep(1L,k3)),nrow=n3,ncol=ncol(ans$g)) else matrix(sample(0L:1L,n3*ncol(ans$g),replace=TRUE,prob=c(1-r,r)),nrow=n3)
                      })
    }
    return(ans)
}

## n1=3
## n=5
## B=10
## expect_warning(s1 <- simulate_trials(B,n1,n))
## expect_equal(length(s1$x),n)
## expect_equal(nrow(s1$g),n)
## expect_equal(ncol(s1$g),pmin(prod(choose(c(n1,n-n1),floor(1/2*c(n1,n-n1)))),B))
## expect_equal(unique(colSums(s1$g)),sum(floor(1/2*c(n1,n-n1))))
## expect_equal(unique(colSums(s1$g[1:n1,])),floor(1/2*n1))
## expect_equal(unique(colSums(s1$g[(n1+1):n,])),floor(1/2*(n-n1)))

##' @title sample size formula for z-test (with known variance)
##' @param delta mean
##' @param sd standard deviation
##' @param sig.level significance level
##' @param power desired power
##' @param type either \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}
##' @param alternative either \code{"one.sided"} or \code{"two.sided"}
##' @return sample size (rounded to the next larger integer)
##' @author Florian Klinglmueller
##' @examples
##' power.z.test(delta=1,sd=1,power=.8)
##' power.t.test(delta=1,sd=1,power=.8)
##' @export
power.z.test <- function(delta,sd,sig.level=.025,power=0.8,
                         type = c("two.sample", "one.sample", "paired"),
                         alternative = c("one.sided", "two.sided")){
    type=match.arg(type)[1]
    alternative=match.arg(alternative)[1]
    if(alternative == "two.sided") sig.level  <- sig.level/2
    dfact <- qnorm(sig.level,lower.tail=F)+qnorm(1-power,lower.tail=F)
    if(type == "two.sample") dfact <- dfact*sqrt(2)
    cat("NOTE: n is the number in *each* group\n")
    ceiling((sd/delta)^2 * (dfact)^2)
}

##' @title sample size formula for the (one-sided) Mann-Whitney U test
##' @param res relative effect size, i.e. probability that Y_i > X_i
##' @param propn proportion of control group samples among total sample
##' @param sig.level significance level
##' @param power desired power
##' @param silent should explanations be suppressed
##' @param n total sample size
##' @return sample size (rounded to the next larger integer)
##' @author Florian Klinglmueller
##' @examples
##' power.z.test(delta=1,sd=1,power=.84)
##' power.t.test(delta=1,sd=1,power=.84)
##' power.u.test(delta2res(delta=1,sd=1),power=.8)
##' @references
##' Noether, Gottfried E. "Sample size determination for some common nonparametric tests." Journal of the American Statistical Association 82.398 (1987): 645-647.
##' @export
power.u.test <- function(res,propn=1/2,sig.level=.025,power=0.8,silent=F){
    N = (qnorm(sig.level,lower.tail=F)+qnorm(1-power,lower.tail=F))^2/(12*propn*(1-propn)*(res-1/2)^2)
    if(!silent) cat("NOTE: n is the number in the *control* group\n")
    ceiling(propn*N)
}

##' Convert difference in means to relative effect size for normally distributed samples
##'
##' @title Effect size to relative effect size
##' @param delta difference in means
##' @param sd common standard deviation
##' @return probability that observations in the treatment sample are larger than in the control sample
##' @examples
##' res <- delta2res(1)
##' power.u.test(res)
##' power.z.test(delta=1,sd=1,power=0.84)
##' power.t.test(delta=1,sd=1,power=0.84)
##' @author float
##' @export
delta2res <- function(delta,sd=1){
    pnorm(0,delta,sqrt(2*sd^2),lower.tail=F)
}

    
##' Conditional power rule for the one-sample (or paired) z-test using the normal distribution sample size formula. Reestimates the standard deviation from the first stage and recomputes the sample size such that the power to reject the null meets the target power assuming that the mean (paired treatment difference) is equal to a prespecified value.
##'
##' @title Conditional power sample size reassessment rule (one-sample z-test)
##' @template power_rules
##' @author Florian Klinglmueller
##' @export
cond_power_rule_norm <- function(x1,m=1,target=.9,alpha=.025,maxN=Inf){
    s <- var(x1)
    dfact <- qnorm(alpha,lower.tail=F)+qnorm(1-target,lower.tail=F)
    min(ceiling(s/m^2 * (dfact)^2),maxN)
}

##' Conditional power rule for the one-sample (or paired) t-test using the function \code{link{stats::power.t.test}}. Reestimates the standard deviation from the first stage and recomputes the sample size such that the power to reject the null meets the target power assuming that the mean (paired treatment difference) is equal to a prespecified value.
##' 
##' @title Conditional power sample size reassessment rule (one-sample t-test)
##' @template power_rules
##' @author Florian Klinglmueller
##' @export
cond_power_rule_t <- function(x1,m=1,target=.9,alpha=.025,maxN=Inf){
    min(maxN,ceiling(power.t.test(power=target,delta=m,sd=sd(x),sig.level=alpha,type='one.sample',alternative='one.sided')$n))
}


##' Compute the robust scale measure Sn given in "Rousseeuw and Croux (1993)" which is about 58% efficient efficient compared to the standard deviation and has a breakdown point of 50%.
##'
##' The default scaling factor is chosen such that the estimate is approximately consistent for normal observations. For different reference distributions other factors may be more apropriate. Note however, that we have hardcoded the bias correction constants applied to small samples - these may differ as well for other reference distributions. Changing the factor may be futile.
##'
##' Note that we use a computationally suboptimal implementation that uses \code{O(n^2)} operations. If time allows we may implement the more efficient method given in "Croux and Rousseeuw (1992)" which requires only \code{O(n*log(n))} operations.
##' @title Robust scale estimate Sn
##' @param x vector of observations
##' @param factor Scaling factor
##' @return numeric value 
##' @author Florian Klinglmueller
##' @export
sn <- function(x,factor=1.1926){
    n <- length(x)
    cns <- c(0.743,1.851,0.954,1.351,0.993,1.198,1.005,1.131)
    cn <- ifelse(n>9,n/(n-(n%%2/10)),cns[n-1])
    factor*sort(matrixStats::rowOrderStats(abs(matrix(x,n,n) - matrix(x,n,n,byrow=T)),which=floor(n/2)+1))[floor((n+1)/2)]
}

##' Compute the robust scale measure Qn given in "Rousseeuw and Croux (1993)" which is about 86% efficient compared to the standard deviation and has a breakdown point of 50%.
##'
##' The default scaling factor is chosen such that the estimate is approximately consistent for normal observations. For different reference distributions other factors may be more apropriate. Note however, that we have hardcoded the bias correction constants applied to small samples - these may differ as well for other reference distributions. Changing the factor may be futile.
##'
##' Note that we use a computationally suboptimal implementation that uses \code{O(n^2)} operations. If time allows we may implement the more efficient method given in "Croux and Rousseeuw (1992)" which requires only \code{O(n*log(n))} operations.
##' @title Robust scale estimate Qn
##' @param x vector of observations
##' @param factor Scaling factor
##' @return numeric value 
##' @author Florian Klinglmueller
##' @export
qn <- function(x,factor=2.2219){
    n <- length(x)
    cns <- c(0.399,0.994,0.512,0.844,0.611,0.857,0.669,0.872)
    cn <- ifelse(n>9,n/(n+(3.8-n%%2*2.4)),cns[n-1])
    xm <- abs(matrix(x,n,n) - matrix(x,n,n,byrow=T))
    factor*sort(xm[upper.tri(xm)])[choose(floor(n/2)+1,2)]
}


##' Compute the robust scale measure Qn pooled over two samples.
##'
##' The default scaling factor is chosen such that the estimate is approximately consistent for normal observations. For different reference distributions other factors may be more apropriate. Note however, that we have hardcoded the bias correction constants applied to small samples - these may differ as well for other reference distributions. Changing the factor may be futile.
##'
##' Note that we use a computationally suboptimal implementation that uses \code{O(n^2)} operations. If time allows we may implement the more efficient method given in "Croux and Rousseeuw (1992)" which requires only \code{O(n*log(n))} operations.
##' @title Robust scale estimate Qn pooled sample
##' @param x vector of observations
##' @param y vector of observations
##' @param factor Scaling factor
##' @return numeric value 
##' @author Florian Klinglmueller
##' @export
qnp <- function(x,y,factor=2.2219){
    n1 <- length(x)
    n2 <- length(y)
    cns <- c(0.249,0.699,0.451,0.773,0.577,0.823,0.657,0.856,0.712,0.879)
    cn <- ifelse(n1>11,n1/(n1 + (n1%%2) * 1.515 + (1-n1%%2) * 3.883),cns[n1-1])
    if(n1 != n2) warning('Bias correction factor may not be correct for unequal sample sizes')
    xm <- abs(matrix(x,n1,n1) - matrix(x,n1,n1,byrow=T))
    ym <- abs(matrix(y,n2,n2) - matrix(y,n2,n2,byrow=T))
    cn*factor*sort(c(xm[upper.tri(xm)],ym[upper.tri(ym)]))[choose(floor(n1/2)+1,2)+choose(floor(n2/2)+1,2)]
}


##' Compute the robust scale measure Sn pooled over two samples. 
##'
##' The default scaling factor is chosen such that the estimate is approximately consistent for normal observations. For different reference distributions other factors may be more apropriate. Note however, that we have hardcoded the bias correction constants applied to small samples - these may differ as well for other reference distributions. Changing the factor may be futile.
##'
##' Note that we use a computationally suboptimal implementation that uses \code{O(n^2)} operations. If time allows we may implement the more efficient method given in "Croux and Rousseeuw (1992)" which requires only \code{O(n*log(n))} operations.
##' @title Robust scale estimate Qn pooled sample
##' @param x vector of observations
##' @param y vector of observations
##' @param factor Scaling factor
##' @return numeric value 
##' @author Florian Klinglmueller
##' @export
snp <- function(x,y,factor=1.1926){
    n1 <- length(x)
    n2 <- length(y)
    cns <- c(0.982,1.302,0.846,1.138,0.877,1.078,0.899,1.048,0.914,1.032)
    cn <- ifelse(n1>11,n1/(n1 + (1-n1%%2) * 1.0326 - (n1%%2) * 0.2469),cns[n1-1])
    if(n1 != n2) stop('Unequal sample sizes not yet implemented')
    xm <- matrixStats::rowOrderStats(abs(matrix(x,n1,n1) - matrix(x,n1,n1,byrow=T)),which=floor(n1/2)+1)
    ym <- matrixStats::rowOrderStats(abs(matrix(y,n2,n2) - matrix(y,n2,n2,byrow=T)),which=floor(n2/2)+1)
    cn*factor*sort(c(xm,ym))[n1]
}

##' Compute the median absolute deviation pooled over two samples
##'
##' We subtract the median from each group and then compute the MAD from the combination of the (shifted) samples. Finite sample corretion factors have been estimated from a simulation study and work well to produce nearly unbiased estimates for normally distributed samples.
##' 
##' The scaling factor does nothing.
##'
##' @title MAD of pooled samples
##' @param x vector of observations
##' @param y vector of observations
##' @param factor Scaling factor (deprecated)
##' @return numeric value 
##' @author Florian Klinglmueller
##' @export
madp <- function(x,y,factor=1){
    n1 <- length(x)
    cns <- c(1.054,1.429,1.239,1.178,1.167,1.111,1.111,1.080,1.081,1.063)
    cn <- ifelse(n1>11,n1/(n1 - 0.659),cns[n1-1])
    cn*mad(c(x-median(x),y-median(y)))
}


##' Compute the interquartile range pooled over two samples
##'
##' We subtract the median from each group and then compute the interquartile range from the combination of the (shifted) samples. Finite sample corretion factors have been estimated from a simulation study. Those do not work particularly well to produce nearly unbiased estimates for normally distributed samples but are better than nothing.
##'
##' The scaling factor does nothing.
##'
##' @title MAD of pooled samples
##' @param x vector of observations
##' @param y vector of observations
##' @param factor Scaling factor (deprecated)
##' @return numeric value 
##' @author Florian Klinglmueller
##' @export
iqrp <- function(x,y,factor=1){
    n1 <- length(x)
    cns <- c(1.186,1.241,1.197,1.173,1.159,1.136,1.129,1.112,1.108,1.097)
    cn <- ifelse(n1>11,n1/(n1 - 1.046),cns[n1-1])
    cn*IQR(c(x-median(x),y-median(y)))/(2*qnorm(3/4))
}



##' @title Robust pooled variance estimate 
##' @param x control group observations
##' @param y treatment group observations
##' @param type what estimator to use
##' @param factor multiplication constant
##' @return robust variance estimate
##' @author float
##' @export
robust_pooled_variance <- function(x,y,type=c('qn','sn','iqr','mad'),factor=NULL){
    type  <- match.arg(type)
    scale_m <- switch(type,
                      'qn' = function(xs,ys) qnp(x,y,factor=ifelse(is.null(factor),2.2219,factor)),
                      'sn' = function(xs,ys) snp(x,y,factor=ifelse(is.null(factor),1.1926,factor)),
                      'iqr' = function(xs,ys) iqrp(x,y)*ifelse(is.null(factor),1,factor),
                      'mad' = function(xs,ys) madp(x,y,factor=ifelse(is.null(factor),1,factor)))
    scale_m(x,y)^2
}


##' Conditional power rule for the two-sample t-test using the function using the normal distribution sample size formula. Reestimates the standard deviation from the first stage and recomputes the sample size such that the power to reject the null meets the target power assuming that the mean (paired treatment difference) is equal to a prespecified value.
##' 
##' @title Conditional power sample size reassessment rule (two-sample z-test)
##' @template power_rules_ts
##' @param ... additional arguments to \code{\link{robust_pooled_variance}}
##' @author Florian Klinglmueller
##' @export
cond_power_rule_norm_ts <- function(x1,y1,delta=1,target=.9,alpha=0.025,maxN=length(x1)*6,rob_var=T,...){
    var <- ifelse(rob_var,
                  robust_pooled_variance(x1,y1,...),
                  pooled_variance(c(x1,y1),rep(0:1,c(length(x1),length(y1)))))
    nE <- 2*(qnorm(alpha,lower=F)+ qnorm(target))^2*var/(delta^2)
    ceiling(min(maxN,nE))
}

##' Conditional power rule for the two-sample t-test using the function using the normal distribution sample size formula. Reestimates the standard deviation from the first stage and recomputes the sample size such that the power to reject the null meets the target power assuming that the mean (paired treatment difference) is equal to a prespecified value.
##' 
##' @title Conditional power sample size reassessment rule (two-sample z-test)
##' @template power_rules_ts
##' @param ... additional arguments to \code{\link{robust_pooled_variance}}
##' @author Florian Klinglmueller
##' @export
cond_power_rule_t_ts <- function(x1,y1,delta=1,target=.9,alpha=0.025,maxN=length(x1)*6,rob_var=T,...){
    var <- ifelse(rob_var,
                  robust_pooled_variance(x1,y1,...),
                  pooled_variance(c(x1,y1),c(rep(0,length(x1)),rep(1,length(y1)))))
    nE <- try(power.t.test(power=target,delta=delta,sd=sqrt(var),sig.level=alpha,type='two.sample',alternative='one.sided')$n,silent=T)
    if(class(nE)=='try-error'){
        nE <- 1+2*(qnorm(alpha,lower=F)+ qnorm(target))^2*var/(delta^2)
        warning(paste0('Numerical fail in t-test SSR, use normal approximation: nE=',nE,' n1=',length(x1)))
    }
    ceiling(min(maxN,nE))
}
##' Conditional power rule for the Mann-Whitney U test, where we estimate the relative effect size based on the observed first stage outcomes by simply counting the number of treatment group samples that have larger outcome values than control group samples. 
##'
##' @title Conditional power sample size reassessment rule (Mann-Whitney U Test)
##' @param x1 First stage control group observations 
##' @param y1 First stage treatment group observations
##' @param target Target power
##' @param alpha Significance level
##' @param maxN Maximum sample size
##' @return Total sample size (in the control group) required to achieve the target power
##' @author float
##' @export
cond_power_rule_u_ts <- function(x1,y1,target=.9,alpha=0.025,maxN=length(x1)*6){
    res <- sum(sapply(y1,function(y) sum(x1 < y)))/(length(x1)*length(y1))
    nE <- power.u.test(res,length(x1)/(length(x1)+length(y1)),alpha,target,silent=T)
    ceiling(min(maxN,nE))
}

##' Computes the inverse normal combination (sqrt(w1)*qnorm(1-p1) + sqrt(w2)*qnorm(1-p2)) of two (independent) p-values
##'
##' @title Inverse normal combination function
##' @param p1 First stage p-value
##' @param p2 Second stage p-value
##' @param w1 First stage weight
##' @param w2 Second stage weight
##' @return p-value corresponding to the inverse normal combination z-score
##' @author Florian Klinglmueller
inverse_normal <- function(p1,p2,w1,w2){
    (sqrt(w1) * qnorm(p1,lower=F) + sqrt(w2) * qnorm(p2,lower=F)) %>% pnorm(lower=FALSE)
}


##' recycle vector to match the length of another
##'
##' @title recycle
##' @param a vector to be recycled
##' @param b vector to be matched
##' @return vector of length at least \code{b}
##' @author Florian Klinglmueller
recycle <- function(a,b){
    if(length(a) < length(b)){
        rep(a,length.out(b))
    } else {
        a
    }
}

##' Implements basic bisection search for root finding. In contrast to uniroot this can also deal with the case where the result of \code{fun} is binary - in this case the \code{bisect} finds the changepoint. 
##'
##' @title Bisection search
##' @param lower lower bound of search range
##' @param upper upper bound of search range
##' @param fun function whose root is to be found
##' @param tol error tolerance
##' @param ... additional arguments to be passed to the function
##' @return root
##' @author float
##' @examples
##' set.seed(5449219)
##' x <- rnorm_cont(16,25,8,cshift=0,csd=24)
##' y <- rnorm_cont(16,19,8,cshift=0,csd=24)
##' cond_power_rule_norm_ts(x[1:8],y[1:8],delta=6,target=.8,rob_var=F)
##' nE <- cond_power_rule_norm_ts(x[1:8],y[1:8],delta=6,target=.8,rob_var=T)
##' xE <- rnorm_cont(nE-16,25,8,cshift=0,csd=24)
##' yE <- rnorm_cont(nE-16,19,8,cshift=0,csd=24)
##' 
##'
##' combtest <- function(delta,alpha){
##'    adaptive_invnormtest_2s(c(y,yE)+delta,c(x,xE),8,16,nE,alpha=alpha)
##' }
##' 
##' bisect(-10, 10, combtest, alpha = 0.5)
##' @export
bisect <- function(lower,upper,fun,tol=0.0001,...){
    opts <- list(...)
    mid  <- (upper + lower)/2
    if(upper < lower){
        tmp <- lower
        lower <- upper
        upper <- tmp
        warning("Lower bound larger than upper bound")
    }
    if(abs(upper - lower) < tol){
        return(mid)
    }
    lambda <- function(lower,upper,fun,flow,fup,opts){
        mid  <- {upper + lower}/2
        if(upper - lower < tol){
            return(mid)
        }
        fmid <- do.call(fun,c(mid,opts))
        if(sign(fmid) == sign(flow)){
            flow  <- fmid
            lower  <- mid
        } else {
            fup  <- fmid
            upper  <- mid
        }
        lambda(lower,upper,fun,flow,fup,opts)
    }
    flow <- do.call(fun,c(lower,opts))
    fup <- do.call(fun,c(upper,opts))
    if(sign(flow) == sign(fup)){
        stop('not different')
    }
    lambda(lower,upper,fun,flow,fup,opts)
}

