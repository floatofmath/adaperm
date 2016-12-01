##' This function computes the three probabilities P(X < Y), P(X < Y and X < Y'), and P(X < Y and X' < Y) [see (2.16), (2.19) and (2.20) on p. 69f of Lehmann (2006)] for two given samples of observations.
##'
##' @title Pair-wise order probabilities
##' @param x vector of (control group) observations
##' @param y vector of (treatment group) observations
##' @return list of order probabilities
##' @export order_probabilities
##' @author float
order_probabilities <- function(x,y){
    m <- length(x)
    n <- length(y)
    mn <- m*n
    p1 <- sum(outer(x,y,FUN='<'))/mn
    o <- outer(y,y,FUN=pmin)
    diag(o) <- -Inf
    p2 <- sum(outer(x,o,FUN='<'))/(mn*(n-1))
    o <- outer(x,x,FUN=pmax)
    diag(o) <- Inf
    p3 <- sum(outer(y,o,FUN='>'))/(mn*(m-1))
    return(list('p1'=p1,'p2'=p2,'p3'=p3))
}


##' This function computes the power and sample size according to the asymptotic power formula given in [(2.23) p.71 of Lehmann (2006)]. If \code{is.null(n)} the sample size required to achieve the target power is computed using bisection search.
##'
##' @title Power and sample size of the Wilcoxon Mann Whitney U test
##' @param n Sample size in the control group
##' @param p1 probability that X < Y
##' @param p2 probability that X < Y and X < Y'
##' @param p3 probability that X < Y and X' < Y
##' @param propn proportion of the total number of observations in the control group
##' @param sig.level significance level of the test
##' @param power desired target power
##' @param silent should hints be suppressed
##' @return Either the required number of control group observations to achieve the targe power, or the power of the WMW test under the specified alternative and given sample size.
##' @export power.w.test
##' @author float
power.w.test <- function(n=NULL,p1=NULL,p2=NULL,p3=NULL,propn=1/2,sig.level=0.025,power = 0.9,silent=F){
    if(length(p1) == 3) {
        p2 <- p1[[2]]
        p3 <- p1[[3]]
        p1  <- p1[[1]]
    }
    pw <- function(m,p1,p2,p3,alpha,r){
        n <- m/r - m
        ## if(any(n*m*{p1*(1-p1) + (n-1)*(p2-p1^2) + (m-1)*(p3-p1^2)} < 0)) {
        ##     warning(paste('Variance estimate not positive:',p1,p2,p3))
        ##     return(ifelse(p1>1/2,1,0))
        ## }
        pnorm({sqrt(n*m*(n+m+1)/12) * qnorm(alpha,lower=F) + (n*m-1)/2 - n*m*p1}/
              sqrt(n*m*{p1*(1-p1) + (n-1)*(p2-p1^2) + (m-1)*(p3-p1^2)}),lower=F)
    }
    if(is.null(n)){
        ceiling(bisect(2,20000,function(m) pw(m,p1,p2,p3,alpha=sig.level,propn) >= power,tol=.3))
    } else {
        pw(n,p1,p2,p3,sig.level,propn)
    }
}

##' Sample size reassessment rule for rank based randomization tests
##'
##' Re-computes the sample size based on first stage observations and assuming that a WMW test is performed at the end of the trial. It estimates the overall (across all stages) number of control-group observations required to achieve the target power. The returned value may be smaller than the preplanned sample size, even the first stage sample size - so any minimum sample size restrictions need to be enforced outside of the function. The maximum (overall control-group) sample size however is enforced internally.
##'
##' The function \code{cond_power_rule_w_ts} recomputes the sample size using the sample size formula of \code{\link{power.w.test}} to achieve the target power \code{target} assuming the order probabilities \code{pp} with a second stage WMW test at the level of the conditional error rate of the preplanned test.
##'
##' The function \code{predictive_power_rule_w_ts} recomputes the order probabilities based on the first stage observations and reestimates the sample size required to achieve the target power using the sample size formula of \code{\link{power.w.test}}. 
##'
##' The functions \code{approximate_power_rule_w_ts}, \code{combination_power_rule_w_ts} and \code{optimal_power_rule_w_ts} compute the sample size that maximizes the combined objective \code{CP(x1,y1) - lambda*nA} defined in [Jennison and Turnbull (2015)]. \code{approximate_power_rule_w} estimates the conditional power using the power formula of \code{\link{power.w.test}} with weighted averages of the prespecified order proabibilities \code{pp} and the order probabilities estimated from the first stage obseravtions. \code{comb_power_rule_w_ts} estimates the conditional power using the power formula of \code{\link{power.w.test}} of a second stage WMW with the prespecified order probabilities at the level of the conditional error of the inverse normal combination WMW test, finally \link{optimal_power_rule_w_ts} estimates the conditional power in the same way however using the conditional error rate of the preplanned randomization test - which requires knowledge of the blinded second stage observations.
##'
##' @name SSR_WMW
##' @param x1 first stage (control-group) sample 
##' @param y1 first stage (control-group) sample
##' @param z2 second stage blinded combined sample
##' @param pp list with prespecified order probabilities p1, p2, p3 (see \code{\link{order_probabilities}} and \code{\link{delta2res}} for details)
##' @param lambda penalty factor for additional sample in the combined objective
##' @param alpha significance level of the preplanned test 
##' @param maxN maximum overall (control-group) sample size
##' @param propn proportion of the total number of observations in the control group
##' @param m preplanned overall (control-goup) sample size
##' @param target desired target power
##' @seealso \code{\link{order_probabilities}} that estimates order probabilities from observations, \code{\link{delta2res}} that computes order probabilisities for normally distributed observations
##' @references Lehmann, Erich Leo, and H. J. D'abrera. Nonparametrics: statistical methods based on ranks. Springer, 2006.
##'
##'             Jennison, Christopher, and Bruce W. Turnbull. "Adaptive sample size modification in clinical trials: start small then ask for more?." Statistics in medicine 34.29 (2015): 3793-3810.
##'
##'             Bauer, Peter, and Franz Koenig. "The reassessment of trial perspectives from interim data—a critical view." Statistics in medicine 25.1 (2006): 23-36.
##' @return reassessed overall sample size
##' @author float
NULL
#> NULL


##' @describeIn SSR_WMW Conditional power rule for the WMW test
##' @export cond_power_rule_w_ts
cond_power_rule_w_ts <- function(x1,y1,z2,pp,target = 0.9,alpha=0.025,maxN=length(x1)*6,propn=1/2){
    m1 <- length(x1)
    m2 <- length(z2)*propn
    n2 <- length(z2) - m2
    A <- permutation_cer(c(x1,y1),z2,rep(0:1,c(length(x1),length(y1))),nt2=n2,test_statistic=ranksum,permutations=10000,alpha=alpha,restricted=TRUE,cer_type='randomized',stratified=TRUE)
    if(A == 0) return(m1+m2)
    if(power.w.test(n=2,p1=pp,sig.level=A,propn=propn,silent=T) >= target) return(m1+2)
    ## psc <- order_probabilities(x1,y1)
    min(m1+power.w.test(p1=pp,sig.level=A,power=target,silent=T,propn=propn),maxN)
}



##' @describeIn SSR_WMW Predicitive power rule for the WMW test
##' @export predictive_power_rule_w_ts
predictive_power_rule_w_ts <- function(x1,y1,m=2*length(x1),target = 0.9,alpha=0.025,maxN=length(x1)*6,propn=1/2){
    m1 <- length(x1)
    n1 <- length(y1)
    mn1 <- length(x1)*length(y1)
    ps <- order_probabilities(x1,y1)
    for(i in names(ps)) assign(i,ps[[i]])
    pw <- function(m,p1,p2,p3,alpha,r){
        n <- m/r-m
        if(any(n*m*{p1*(1-p1) + (n-1)*(p2-p1^2) + (m-1)*(p3-p1^2)} < 0)) {
            warning(paste('Variance estimate not positive:',p1,p2,p3))
            return(ifelse(p1>1/2,1,0))
        }
        pnorm({sqrt(n*m*(n+m+1)/12) * qnorm(alpha,lower=F) + (n*m-1)/2 - n*m*p1}/
              sqrt(n*m*{p1*(1-p1) + (n-1)*(p2-p1^2) + (m-1)*(p3-p1^2)}),lower=F)
    }
    if(power.w.test(m1,p1=ps,sig.level=alpha,propn=propn)>target) return(m1)
    if(power.w.test(maxN,p1=ps,sig.level=alpha,propn=propn)<target) return(maxN)
    power.w.test(p1=ps,sig.level=alpha,propn=propn,power=target)
}


##' @describeIn SSR_WMW Efficient power rule for the WMW test using approximate conditional order probabilities
##' @export approximate_power_rule_w_ts
approximate_power_rule_w_ts <- function(x1,y1,m=2*length(x1),pp,lambda = 1e-4,alpha=0.025,maxN=length(x1)*6,propn=1/2){
    m1 <- length(x1)
    m2 <- m - m1
    ps <- order_probabilities(x1,y1)
    nA <- 0:(maxN-(m1+m2))
    pc <- lapply(1:3,function(i) {(m1*ps[[i]] + (m2+nA)*pp[[i]])/(m1+m2+nA)})
    cp <- power.w.test(m2+nA,p1=pc,sig.level=alpha,propn=propn) - lambda * 0:(maxN-(m1+m2))
    mA <- which.max(cp)-1
    return(c(m1+m2+mA))
}

##' @describeIn SSR_WMW Efficient power rule for the WMW test using the combination test conditional power
##' @export combination_power_rule_w_ts
combination_power_rule_w_ts <- function(x1,y1,m=2*length(x1),pp,lambda = 1e-4,alpha=0.025,maxN=length(x1)*6,propn=1/2){
    m1 <- length(x1)
    m2 <- m - m1
    A <- pnorm({qnorm(alpha,lower=F) - sqrt(m1/m)*qnorm(wilcox.test(y1,x1,alternative='greater')$p.value,lower=F)}/sqrt(1-m1/m),lower=F)
    cp <- power.w.test(m2+0:(maxN-(m1+m2)),p1=pp,sig.level=A,propn=propn) - lambda * 0:(maxN-(m1+m2))
    mA <- which.max(cp)-1
    return(c(m1+m2+mA))
}

##' @describeIn SSR_WMW Efficient power rule for the WMW test using the combination test conditional power
##' @export combination_power_rule_w_ts
combination_power_rule_w_ts_f <- function(x1,y1,m=2*length(x1),pp,lambda = 1e-4,alpha=0.025,maxN=length(x1)*6,propn=1/2,futility=0.8){
    m1 <- length(x1)
    m2 <- m - m1
    A <- pnorm({qnorm(alpha,lower=F) - sqrt(m1/m)*qnorm(wilcox.test(y1,x1,alternative='greater')$p.value,lower=F)}/sqrt(1-m1/m),lower=F)
    cp <- power.w.test(m2+0:(maxN-(m1+m2)),p1=pp,sig.level=A,propn=propn) - lambda * 0:(maxN-(m1+m2))
    mA <- which.max(cp)-1
    ## note to self any ideas of using lambda to control the futility boundary do not appear very practical. a) negative futility at preplanned sn does not work as cp - l*0 >= 0 b) changing to l*nE-n1 does not work either because for the studied cases the maximum of cp - l*n is always somewhere n<m2
    if(max(cp+lambda * 0:(maxN - (m1 + m2)))<=futility) return(NULL)
    return(c(m1+m2+mA))
}

##' @describeIn SSR_WMW Efficient power rule for the WMW test using the conditional power of a second stage WMW test
##' @export optimal_power_rule_w_ts
optimal_power_rule_w_ts <- function(x1,y1,z2,pp,lambda = 1e-4,alpha=0.025,maxN=length(x1)*6,propn=1/2){
    m1 <- length(x1)
    m2 <- length(z2)*propn
    n2 <- length(z2) - m2
    A <- permutation_cer(c(x1,y1),z2,rep(0:1,c(length(x1),length(y1))),nt2=n2,test_statistic=ranksum,permutations=10000,alpha=0.025,restricted=TRUE,cer_type='randomized',stratified=TRUE)
    nA <- 0:(maxN-(m1+m2))
    cp <- power.w.test(m2+nA,p1=pp,sig.level=A,propn=propn) - lambda*nA
    mA <- which.max(cp)-1
    return(c(m1+m2+mA))
}

lambda <- list()
lambda[[5]] <- 345e-4
lambda[[15]] <- 125e-4
lambda[[50]] <- 395e-5

lambdas <- c(345e-4,125e-4,395e-5,450e-4,160e-4,550e-5)
sigmas <- c(0.85,1.55,2.8,0.8,1.5,2.6)
n1s <- c(5,15,50,5,15,50)
lambda.model <- lm(log(lambdas)~sigmas+n1s)
lambdas.fit <- exp(-1.8 - 1.9 * sigmas + 0.03 * n1s)

## exp(predict(lambda.model,newdata=list(sigmas=1.5,n1s=15)))
## plot(lambdas,lambdas.fit,xlim=c(0,.05),ylim=c(0,.05))
## abline(0,1)

##' Estimate the penalty factor lambda from n1, in a way to give roughly 80% overall power
##'
##' Except for first stage sample sizes 5, 15 and 50 this may give totally useless results.
##'
##' @title n1 to lambda
##' @param n1 first stage control group sample size
##' @param target target power
##' @return penalty factor for the combined sample size objective
##' @seealso \code{\link{optimal_power_rule_w_ts}}, \code{\link{combination_power_rule_w_ts}} and  \code{\link{approximate_power_rule_w_ts}} for sample size rules based on a combined objective (maximum conditional power with moderate small sample size)
##' @export
##' @author float
n2lambda <- function(n1,sigma,target=.8){
    if(target != .8) stop('Any other target than 80% power is not yet supported')
    if(!(n1 %in% n1s && sigma %in% sigmas)){
        return(predict(lambda.model,newdata=list(n1s=n1,sigmas=sigma)))
    }
    lambdas[which(n1s %in% n1&sigmas %in% sigma)]
}

    
