##' Quick no-nonsens adaptive permutation decision rule
##' 
##' @title Quick adaptive permutation test
##' @param x1 first stage observations
##' @param x2 second stage observations
##' @param xE extended stage observations
##' @param g1 frist stage treatment assignments
##' @param g2 second stage treatment assignments
##' @param gE extended stage treatment assignments
##' @param test_statistic test statistic to use 
##' @param alpha significance level
##' @param permutations number of permutations to use
##' @param restricted should the treatment group sizes be fixed
##' @param atest_type if 'CER' compute only conditional error rate, else type of adaptive test should be performed (see \code{\link{perm_test}} for details) 
##' @param cer_type what type of conditional error rate function should be used (see \code{\link{permutation_cer}}) for details
##' @return pvalue
##' @author Florian Klinglmueller
##' @export
adaptive_permdr <- function(x1,x2,xE,
                        g1,g2,gE,
                        test_statistic,
                        alpha,
                        permutations,
                        restricted,
                        atest_type=c("non-randomized","mid-p","davison_hinkley","CER"),
                            cer_type=c("non-randomized","randomized","uniform")){
    if(length(xE)==0 && !atest_type=='CER'){
        return(alpha>=perm_test(x1,x2,g1,g2,stat=test_statistic,B=permutations,restricted=restricted,type=atest_type))
    }
    A <- permutation_cer(x1,x2,
                         g1,sum(g2>0),
                         test_statistic=test_statistic,
                         alpha=alpha,
                         permutations=permutations,
                         restricted=restricted,
                         cer_type=cer_type)
    if(atest_type == 'CER'){
        return(A)
    } else {
        return(A >= perm_test(x2,xE,g2,gE,stat=test_statistic,B=permutations,restricted=restricted,type=atest_type))
    }
}
.cer_types = c("non-randomized","randomized","uniform")
.atest_types = c("non-randomized","mid-p","davison_hinkley","CER")


##' User friendly wrapper to \code{\link{adaptive_permdr}}
##'
##' @title Permutation based decision rule for adaptive designs
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
##' @return Decision \code{TRUE} if null hypothesis is rejected
##' @author Florian Klinglmueller
##' @export
adaperm_DR <- function(x,g=NULL,n1,n,m1=n1,m=n,test_statistic,alpha=.025,cer_type=.cer_types,atest_type=.atest_types,permutations=10000){
    if(is.null(g)){
        ## one-sample test
        restricted <- FALSE
        obs <- split_sample_os(x,n1,n)
    } else {
        restricted  <- TRUE
        if(length(g) == m1+n1 && cer_type=='CER'){
            ## if we only want to know the CER
            g <- c(g,rep(0,n-n1),rep(1,m-m1))
        }        
        obs <- split_sample_ts(x[g<=0],x[g>0],n1,n,m1,m)
    }
    xs <- obs[[1]]
    gs <- obs[[2]]
    adaptive_permdr(xs[[1]],xs[[2]],xs[[3]],
                      gs[[1]],gs[[2]],gs[[3]],
                      test_statistic=test_statistic,
                      alpha=alpha,
                      permutations=permutations,
                      restricted=restricted,
                      atest_type=atest_type,
                      cer_type=cer_type)
}


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



