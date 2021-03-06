
##' Returns a logical matrix with all possible assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control) given that \code{k} objects are in one group
##'
##' @title All combinations of two group assignments
##' @param n number of observations
##' @param k number of observations in (e.g. treatment) group
##' @return integer matrix of size \code{n} x \code{choose(n,k)}
##' @author Florian Klinglmueller
all_reassignments <- function(n,k){
    all_reassignments_cpp(n,k)
    ## N <- choose(n,k)
    ## M <- matrix(0L,n,N)
    ## M[cbind(as.integer(all_combinations_cpp(n,k)),rep(1:N,each=k))] <- 1L
    ## M
}

##' Returns a logical matrix with all possible assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control)
##'
##' @title All assignments
##' @param n number of observations
##' @return integer matrix of size \code{n} x \code{2^n}
##' @author Florian Klinglmueller
all_assignments <- function(n){
    t(e1071::bincombinations(n))
}

##' Returns a logical matrix with B random assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control)
##'
##' @title Random combinations of two group assignments
##' @param n number of observations
##' @param B number of random combinations
##' @return integer matrix of size \code{n} x \code{B}
##' @author Florian Klinglmueller
random_assignments <- function(n,B){
    matrix(sample(0L:1L,n*B,rep=T),n,B)
}

##' Returns a logical matrix with B random assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control) given that \code{k} objects are in one group
##'
##' @title Random reassignments to two groups
##' @param n number of observations
##' @param k number of observations in (e.g. treatment) group
##' @param B number of random combinations
##' @return integer matrix of size \code{n} x \code{B}
##' @author Florian Klinglmueller
random_reassignments <- function(n,k,B){
    random_reassignments_cpp(c(rep(0L,n-k),rep(1L,k)),B)
}
##library(microbenchmark)

## only faster for small c1
## join_to_row1 <- function(c1,c2){
##     do.call('cbind',lapply(1:ncol(c1),function(i) rbind(matrix(c1[,i],nc=ncol(c2),nr=nrow(c1)),c2)))
## }


## slightly slower maybe due to magrittr
## join_to_row2 <- function(c1,c2){
##     list(x=c1,y=c2) %>%
##         lapply(FUN=function(combs) { combs %>% t  %>% as.data.frame %>% mutate(idx=1)}) %$%
##             full_join(x,y,by='idx') %>% mutate(idx=NULL)
## }

##' Joins two matrices repeats each row of the first matrix by the
##' number of rows of the second matrix and column binds them together to a 
##' 
##' @title Join two matrices by row
##' @param c1 a n1 \times m1 matrix
##' @param c2 a n2 \times m2 matrix
##' @return n1 + n2 \times m1*m2 matrix
##' @author Florian Klinglmueller
join_to_row <- function(c1,c2){
    out <- dplyr::mutate(dplyr::full_join(dplyr::mutate(as.data.frame(c1),idx=1),
                                          dplyr::mutate(as.data.frame(c2),idx=1),by='idx'),idx=NULL)
    dimnames(out) <- list(1:nrow(out),1:ncol(out))
    as.matrix(out)
}




## c1 <- resamplingMCP:::all_reassignments(10,7)
## c2 <- resamplingMCP:::all_reassignments(10,5)
## microbenchmark(join_to_row(c1,c2),join_to_row1(c1,c2))

##' Returns all possible assignments to two groups stratified by stages.
##'
##' Number of observations \code{ns} and observations from one group \code{ks} will be recycled to match the longer argument. 
##'
##' @title Stagewise two-group assignments
##' @param ns vector of number of observations per stage
##' @param ks vector of number of observations in group (e.g. treatment) per stage
##' @param restricted 
##' @param B number of assignments
##' @return logical matrix of dimension \code{sum(ns)} \times \code{prod(choose(ns,ks))}
##' @author Florian Klinglmueller
strat_reassignments <- function(ns,ks,restricted=TRUE,B=NULL){
    ns <- recycle(ns,ks)
    ks <- recycle(ks,ns)
    if(restricted){
        n_combs <- prod(choose(ns,ks))
    } else {
        n_combs <- prod(2^ns)
    }
    random <- ifelse(is.null(B),FALSE,(n_combs > B))
    if(random) {
        combinations <- if(restricted) {
            function(n,k) random_reassignments(n,k,B)
        } else {
            function(n,k) random_assignments(n,B)
        }
        cs <- lapply(which(ns>0),function(i) combinations(ns[i],ks[i]))
        out <- do.call('rbind',cs)
        out
    } else {
        combinations <- if(restricted) {
            all_reassignments
        } else {
            function(n,k) all_assignments(n)
        }
        cs <- lapply(which(ns>0),function(i) t(combinations(ns[i],ks[i])))
        t(Reduce(join_to_row,cs))
    }
}

##' Reference/permutation space of a two-stage adaptive permutation test
##'
##' @title Reference space
##' @param g1 First stage treatment assignments
##' @param g2 Second stage treatment assignments
##' @param g3 Third stage treatment assignments
##' @param restricted Whether group sizes are considered fixed
##' @param B Number of permutations to be used (if smaller than all permutations)
##' @param add_obs should the identity permutation be added
##' @return integer matrix 
##' @author Florian Klinglmueller
##'
##' @export
omega <- function(g1,g2=NULL,g3=NULL,restricted = TRUE,B=1000,add_obs=TRUE){
    ns <- sapply(list(g1,g2,g3),length)
    ks  <- sapply(list(g1,g2,g3),function(g) sum(g>0))
    ps <- strat_reassignments(ns,ks,restricted=restricted,B=B)
    n_combs <- ifelse(restricted,prod(choose(ns,ks)),prod(2^ns))
    if(n_combs > B & add_obs) cbind(c(g1,g2,g3),as.matrix(ps)) else as.matrix(ps)
}

##' Compute the conditional permutation distribution given first stage group assignments
##'
##' Second and third stage treatment assignments are only passed to define the second stage sample size and the sizes of the corresponding (treatment) groups. If the number of requested permutations \code{B} is larger than the number of all possible permutations, only the latter will be used.
##'
##' In contrast to \code{\link{perm_dist}} \code{cond_dist} does not add the original test statistic if less than all possible permutations are generated.
##' 
##' @title Conditional permutations distribution
##' @param x1 First stage data
##' @param x2 Second stage data
##' @param g1 First stage treatment assignments
##' @param g2 (Dummy) second stage treatment assignments (see Details)
##' @param stat Function that computes the test statistic (signature \code{function(x,g)})
##' @param B Number of permutations 
##' @param x3 Third stage data (e.g. sample size increase)
##' @param g3 (Dummy) third stage treatment assignments (see Details)
##' @param restricted Should sample sizes be restricted by stage 
##' @param stratified Should permutation be stratified by stage
##' @return numeric vector 
##' @author Florian Klinglmueller
cond_dist <- function(x1,x2,g1,g2,stat,B,x3=NULL,g3=NULL,restricted=TRUE,stratified=TRUE){
    if(stratified){
        omega <- omega(g2,g3,restricted=restricted,B=B,add_obs=FALSE)
    } else {
        omega <- omega(c(g2,g3),restricted=restricted,B=B,add_obs=FALSE)
    }
    omega <- rbind(matrix(g1,nrow=length(g1),ncol=ncol(omega)),omega)
    stat(c(x1,x2,x3),omega)
}



##' Compute the permutation distribution given observations and a test statistic.
##'
##' First, second and third stage treatment assignments are only passed to define the second stage sample size and the sizes of the corresponding (treatment) groups. If the number of permutations is less than the number of all possible permutations the test observed test statistic will be added to the result such that $B+1$ values will be returned. If the number of requested permutations \code{B} is larger than the number of all possible permutations, only the latter will be used. 
##' 
##' @title Permutation distribution
##' @param x1 First stage data
##' @param x2 Second stage data
##' @param g1 (Dummy) first stage treatment assignments
##' @param g2 (Dummy) second stage treatment assignments (see Details)
##' @param stat Function that computes the test statistic (signature \code{function(x,g)})
##' @param B Number of permutations 
##' @param x3 (Dummy) third stage data (e.g. sample size increase)
##' @param g3 (Dummy) third stage treatment assignments (see Details)
##' @param restricted Should group sizes be considered fixed
##' @param stratified Should permutation be stratified by stage
##' @return numeric vector
##' @author Florian Klinglmueller
##' @export
perm_dist <- function(x1,x2,g1,g2,stat,B,x3=NULL,g3=NULL,restricted=TRUE,stratified=TRUE){
    if(stratified){
        omega <- omega(g1,g2,g3,restricted=restricted,B=B)
    } else {
        omega <- omega(c(g1,g2,g3),restricted=restricted,B=B)
    }
    stat(c(x1,x2,x3),omega)
}

##' Permutation test p-value
##'
##' Type "non-randomized" will give the standard p-value, which is define by the proportion of sample permutations with larger or equal test statistic than the observed one. "midp" will give the mid-p value which is the proportion of sample permutations with larger test statistic plus one half the proportion of permutations with equal test statistic. "randomized" will return the proportion of permutations with larger test statistic plus a random number that is uniformly distributed between zero and the proportion of permutations with test statistic equal to the observed.
##' 
##' @title Permutation p-value
##' @param t observed test statistic
##' @param dist permutation distribution
##' @param type type of p-value to compute
##' @return p-value
##' @author Florian Klinglmueller
t2p  <- function(t,dist,type=c('non-randomized','midp','randomized','davison_hinkley')){
    switch(type[1],
           'non-randomized'=mean(dist>=t),
           'midp'=mean(dist>t) + .5*mean(dist==t),
           'randomized'=mean(dist>t) + runif(1)*mean(dist==t),
           'davison_hinkley'=(sum(dist>t)+1)/(sum(dist!=t)+1))
}


##' Level alpha critical boundary of permutation test. Returns the smallest unique value of the permutation distribution such that the proportion of statistics larger is smaller than alpha. 
##' 
##' @title p to t
##' @param p p-value/alpha level
##' @param dist permutation distribution
##' @return critical value
##' @author Florian Klinglmueller
p2t <- function(p,dist){
    max(dist[rank(-dist,ties='max') >= p*length(dist)])
}



##' Perform a permutation test (stratified by stages)
##'
##' \code{type} specifies which permutation p-value to use. "non-randomized" will use the p-value of the non-randomized permutation test. "randomized" will use a randomized p-value that subtracts a random quantity from the non-randomized p-value, such that the a test that rejects for \code{p <= alhpa} has size exactly \code{alpha}. "midp" computes the mid-p-value which is the expected value of the randomized p-value. Davison-Hinkley returns the fraction between the permutation test statistics larger than \code{t} plus one and the permutation test statistics different from \code{t} plus one. The latter is mostly usefull when not all possible permutations are used. 
##' 
##' @title perform permutation test
##' @param x1 First stage observations
##' @param x2 Second stage observations
##' @param g1 First stage group assignments
##' @param g2 Second stage group assignments
##' @param stat Function that computes test statistic 
##' @param B Number of permutations to use
##' @param x3 Third stage observations
##' @param g3 Third stage group assignments
##' @param restricted Should group sizes be considered fixed
##' @param type Type of p-value to compute (see details)
##' @param stratified should permutation be stratified by stage
##' @return p-value of the permutation test
##' @author Florian Klinglmueller
##'
##' @export
perm_test <- function (x1, x2, g1, g2, stat, B, x3 = NULL, g3 = NULL,
                       restricted,
                       type=c('non-randomized','randomized','midp','davison_hinkley'),
                       stratified=TRUE) 
{
    dist <- perm_dist(x1, x2, g1, g2, stat, B, x3, g3, restricted = restricted,stratified=stratified)
    t <- stat(c(x1, x2, x3), c(g1, g2, g3))
    t2p(t,dist,type[1])
}


##' Perform an adaptive permutation test
##'
##' Note that we currently assume that the group sizes of the second
##' stage in the adapted trial are equal to those of the preplanned
##' test.
##' 
##' @title Adaptive permutation test
##' @param x1 First stage observations
##' @param x2 Second stage observations
##' @param g1 First stage treatment assignments 
##' @param g2 Second stage treatment assignments (see Details)
##' @param stat Function that computes the test statistic
##' @param B Number of draws from the permutation distribution to use
##' @param alpha Significance level of the test
##' @param x3 Third stage observations (e.g. sample size increase)
##' @param g3 Third stage treatment allocations (e.g. sample size increase)
##' @return logical indicating wether the test rejects
##'
##' @author Florian Klinglmueller
##' @export
adaptive_perm_test <- function(x1,x2,g1,g2,stat,B,alpha=.025,x3=NULL,g3=NULL){
    cer <- permutation_CER(x1,g1,x2,stat,B,alpha,g2)
    perm_test(x2,x3,g2,g3,stat,B)<=alpha
}

##' Difference between significance level and actual size of the permutation test, due to discreteness
##'
##' @title trimmings
##' @param g1 first stage treatment group assignments
##' @param g2 second stage treatment group assignmens
##' @param restricted should permutation be restricted by group size stratified by stage 
##' @param alpha significance level
##' @return difference between size and significance level
##' @author Florian Klinglmueller
trimmings <- function(g1,g2,restricted,alpha=.025){
    if(!restricted){
        n <- length(g1)+length(g2)
        B <- 2^n
    } else {
        n <- sum(g1<=0)+sum(g2<=0)
        m <- sum(g1>0)+sum(g2>0)
        n1 <- sum(g1<=0)
        m1 <- sum(g1>0)
        B <- choose(m1+n1,m1)*choose(m+n-n1-m1,m-m1)
    }
    b <- floor(B*alpha)
    alpha-(b/B)
}   
