adaptive_permtest_ts <- function(x1,x2,xE,g1,g2,gE,stat,permutations=10000,conf_level=.05){
    pval <- adaptive_permtest_quick(x1,x2,xE,g1,g2,gE,stat,permutations=10000)
    err <- function(d,conf_level=.5){
        x1 <- x1-d*g1
        x2 <- x2-d*g2
        xE <- xE-d*gE
        adaptive_permtest_quick(x1,x2,xE,g1,g2,gE,stat,permutations=10000)-conf_level
    }
    t2 <- stat(c(x2,xE),c(g2,gE))
    names(t2) <- deparse(substitute(stat))
    ## rms
    guesstimate <- sqrt(mean(c(x1,x2,xE)^2))*sign(.5-pval)
    estimate <- c("constant additive effect"=uniroot(err,c(0,guesstimate))$root)
    conf.int <- c(uniroot(err,c(0,guesstimate*sign(conf_level-pval)),conf_level=conf_level)$root,Inf)
    attr(conf.int,"conf.level") <- 1-conf_level
    out <- list(statistic = t2,
                parameter = c('permutations'=permutations),
                p.value = pval,
                estimate = estimate,
                conf.int=conf.int,
                alternative = "one.sided",
                null.value = c("distribution of treated samples"='different from control samples'),
                data.name = paste(deparse(substitute(c(x1,x2,xE,g1,g2,gE))),collapse=', '),
                method = "Two-sample permutation test for adaptive designs")
    class(out) <- 'htest'
    out
}



##' Quick no-nonsens adaptive permutation test that computes only the p-value
##' 
##' @title Quick adaptive permutation test
##' @param x1 first stage observations
##' @param x2 second stage observations
##' @param xE extended stage observations
##' @param g1 frist stage treatment assignments
##' @param g2 second stage treatment assignments
##' @param gE extended stage treatment assignments
##' @param stat test statistic to use 
##' @param permutations number of permutations to use
##' @param restricted should the treatment group sizes be fixed
##' @return pvalue
##' @author Florian Klinglmueller
adaptive_permtest_quick <- function(x1,x2,xE,g1,g2,gE,stat,permutations=10000,restricted=T){
    ## permutation distribution of preplanned test
    pdist  <- perm_dist(x1,x2,g1,g2,stat,B=permutations,restricted=restricted)
    ## conditional distribution of preplanned test given the first stage data
    cdist <- adaperm:::cond_dist(x1,x2,g1,g2,stat,B=permutations,restricted=restricted)
    ## permutation distribution of adapted test
    edist <- perm_dist(x2,xE,g2,gE,stat,B=permutations,restricted=restricted)
    ## observed test statistic of adapted test
    t2 <- stat(c(x2,xE),c(g2,gE))
    ## adaptive p-value 
    pval <- mean(pdist >= quantile(cdist,mean(edist<t2),type=1))
    pval
}


##' Adaptive permutation test one-sample problems
##'
##' @title One-sample adaptive permutation test
##' @template onesample_sims
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
##' @export
adaptive_permtest_os <- function(x,n1,n,ne,test_statistic,perms=50000,alpha=0.025){
    obs <- split_sample_os(x,n1,n,ne)
    if(ne>n){
        xs <- obs[[1]]
        gs <- obs[[2]]
        return(adaptive_permtest_quick(xs[[1]],xs[[2]],xs[[3]],gs[[1]],gs[[2]],gs[[3]],test_statistic,restricted=FALSE,perms))
    } else {
        xs <- obs[[1]]
        gs <- obs[[2]]
        return(perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,restricted=FALSE,B=perms))
    }
}


##' Adaptive decision rule
##'
##' Slightly faster than the test that computes p-values
##' 
##' @title One-sample adaptive permutation test
##' @template onesample_sims
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
adaptive_DR_os <- function(x,n1,n,ne,test_statistic,perms=50000,alpha=0.025){
    obs <- split_sample_os(x,n1,n,ne)
    if(ne>n){
        xs <- obs[[1]]
        gs <- obs[[2]]
        A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,restricted=FALSE,permutations=perms,alpha=alpha)
        q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,restricted=FALSE,B=perms)
        return(A>=q)
    } else {
        xs <- obs[[1]]
        gs <- obs[[2]]
        return(alpha>=perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,restricted=FALSE,B=perms))
    }
}


##' Adaptive decision rule
##'
##' Slightly faster than the test that computes p-values
##' 
##' @title Two-sample adaptive permutation test
##' @template twosample_sims
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
adaptive_DR_ts <- function(x,y,n1,n,ne,m1,m,me,test_statistic,perms=50000,alpha=0.025){
    if(ne>n){
        xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
    } else {
        if(me>m) stop('Stages with controls only not supported')
        xs <- split(x,rep(1:2,c(n1,ne-n1)))
    }
    if(me>m){
        ys <- split(y,rep(1:3,c(m1,m-m1,me-m)))
    } else {
        if(ne>n) stop('Stages with treatments only not supported')
        ys <- split(y,rep(1:2,c(m1,me-m1)))
    }
    gs <- lapply(1:length(xs),function(i) rep(0:1,c(length(xs[[i]],ys[[i]]))))
    xs <- lapply(1:length(xs),function(i) c(xs[[i]],ys[[i]]))
    A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,permutations=perms,alpha=alpha)
    q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,B=perms)
    A>=q
}

##' Adaptive decision rule
##'
##' Slightly faster than the test that computes p-values
##' 
##' @title Two-sample adaptive permutation test
##' @template twosample_sims
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
##' @export
adaptive_permtest_ts <- function(x,y,n1,n,ne,m1=n1,m=m,me=me,test_statistic,perms=50000,alpha=0.025){
    obs <- split_sample_ts(x,y,n1,n,ne,m1,m,me)
    xs <- obs[[1]]
    gs <- obs[[2]]
    if(ne>n){
        return(adaptive_permtest_quick(xs[[1]],xs[[2]],xs[[3]],gs[[1]],gs[[2]],gs[[3]],test_statistic,restricted=TRUE,perms))
    } else {
        xs <- obs[[1]]
        gs <- obs[[2]]
        return(perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,restricted=TRUE,B=perms))
    }
}






##' Adaptive permutation test one-sample problems
##'
##' @title One-sample adaptive permutation test with subsampling
##' @template onesample_sims
##' @param combination_function Function to combine stage-wise (permutation) p-values
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
adaptive_permtest2_os <- function(x,n1,n,ne,test_statistic,perms=50000,resam=100,alpha=0.025){
    obs <- split_sample_os(x,n1,n,ne)
    if(ne>n){
        xs <- obs[[1]]
        gs <- obs[[2]]
        A <- permutation_CER2(xs[[1]],gs[[1]],xs[[2]],xs[[3]],test_statistic,restricted=FALSE,permutations=perms,subsamples=resam,alpha=alpha)
        q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,restricted=FALSE,B=perms)
        return(A>=q)
    } else {
        xs <- obs[[1]]
        gs <- obs[[2]]
        return(alpha>=perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,restricted=FALSE,B=perms))
    }
}



