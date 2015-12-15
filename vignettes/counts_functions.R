
adaptive_permtest_2s <- function(x,y,n1,n,ne,test_statistic,
                                 m1=n1,m=n,me=ne,perms=1000,alpha=0.025){
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
  gs <- lapply(1:length(xs),function(i) rep(0:1,c(length(xs[[i]]),length(ys[[i]]))))
  xs <- lapply(1:length(xs),function(i) c(xs[[i]],ys[[i]]))
  if(me==m){
    return(alpha>=perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,B=perms))
  }
  A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,permutations=perms,alpha=alpha)
  q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,B=perms)
  A>=q
}
###########################

adaptive_permtest2_2s <- function(x,y,n1,n,ne,test_statistic,
                                  m1=n1,m=n,me=ne,perms=1000,alpha=0.025,
                                  resam=100){
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
  gs <- lapply(1:length(xs),function(i) rep(0:1,c(length(xs[[i]]),length(ys[[i]]))))
  xs <- lapply(1:length(xs),function(i) c(xs[[i]],ys[[i]]))
  if(me==m){
    return(alpha>=perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,B=perms))
  }
  A <- permutation_CER2(xs[[1]],gs[[1]],xs[[2]],xs[[3]],test_statistic,one_sample=FALSE,restricted=FALSE,
                        permutations=perms,subsamples=resam,alpha=alpha)
  q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,B=perms)
  A>=q
}

######## other adaptive tests
##################
#' To be checked
adaptive_invnormtest_2s <- function(x,y,n1,n,ne,m1=n1,m=n,me=ne,alpha=0.025){
  xs <- split(x,rep(1:2,c(n1,ne-n1)))
  ys <- split(y,rep(1:2,c(m1,ne-m1)))
  p1 <- t.test(xs[[1]],ys[[1]],alternative='greater')$p.value
  p2 <- t.test(xs[[2]],ys[[2]],alternative='greater')$p.value
  alpha >= {sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)} %>% sum() %>% pnorm(lower=FALSE) 
}

#' To be checked
adaptive_invnormtest_negbin_2s <- function(x,y,n1,n,ne,m1=n1,m=n,me=ne,alpha=0.025){
  xs <- split(x,rep(1:2,c(n1,ne-n1)))
  ys <- split(y,rep(1:2,c(m1,ne-m1)))
  sg1 <- summary(glm.nb(c(xs[[1]],ys[[1]])~rep(0:1,c(n1,m1))))
  sg2 <- summary(glm.nb(c(xs[[2]],ys[[2]])~rep(0:1,c(ne-n1,ne-m1))))
  p1 <- sg1$coefficients[2,"Pr(>|z|)"]
  s1 <- sg1$coefficients[2,"z value"]<0
  p2 <- sg2$coefficients[2,"Pr(>|z|)"]
  s2 <- sg2$coefficients[2,"z value"]<0
  p1 <- ifelse(s1,p1/2,1-p1/2)
  p2 <- ifelse(s2,p2/2,1-p2/2)
  alpha >= {sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)} %>% sum() %>% pnorm(lower=FALSE) 
}

adaptive_waldtest_2s <- function(x,y,n1,n,ne,m1=n1,m=n,me=ne,alpha=0.025){
  log(mean(y)/mean(x))/sqrt(1/sum(x)+1/sum(y))>qnorm(alpha,lower=F)
}
###########


##' @param rule adaptive sample size rule
##' @param rdist random number generator for the data
##' @param test_statistic function that computes the test statistic
##' @param control_opts list of options past to \code{rdist} for the control group
##' @param treatment_opts list of options past to \code{rdist} for the treatment group
##' @param m1 first stage sample size (control group)
##' @param m preplanned total sample size (control group)
##' @param ... 
##' @return list of test results 
##' @author Florian Klinglmueller
compare_adaptive_tests_2s <- function(n1,n,rule,rdist,
                                      test_statistic,
                                      control_opts=list(),
                                      treatment_opts=control_opts,
                                      m1=n1,m=n,resam=100,...){
  x <- do.call(match.fun(rdist),c(n=n,control_opts))
  y <- do.call(match.fun(rdist),c(n=n,treatment_opts))
  nes <- rule(x[1:n1],y[1:m1])
  if(nes[1]>n){
    ne <- nes[1]
    x <- c(x,rdist(ne[1]-n,...))
  } else {
    ne <- n
  }
  if(nes[2]>m){
    me <- nes[2]
    y <- c(x,rdist(ne[2]-m,...))
  } else {
    me <- m
  }
  list(permtest = adaptive_permtest_2s(x,y,n1,n,ne,test_statistic,m1,m,me),
       permtest2 = adaptive_permtest2_2s(x,y,n1,n,ne,test_statistic,m1,m,me,resam=resam),
       invnorm = adaptive_invnormtest_2s(x,y,n1,n,ne,m1,m,me),
       invnorm_nb=adaptive_invnormtest_negbin_2s(x,y,n1,n,ne,m1,m,me),
       wald = adaptive_waldtest_2s(x,y,n1,n,ne))
}

######################################
#OLD (to be removed?)
# run_simulation <- function(n1,n,delta=0,B=10,...){
#     results <- mclapply2(1:B,function(i) 
#       compare_adaptive_tests_2s(n1,n,function(x,y) c(4,4),r_counts,diffmean,treatment_opts=list(delta=delta))
#       )
#     ##:ess-bp-start::conditional@length(unique(sapply(results,length))) != 1:##
# browser(expr={length(unique(sapply(results,length))) != 1})##:ess-bp-end:##
#     results %>% bind_rows %>% colMeans
# }
# 





samplesize_counts <- function(theta,lambda,k=1,t=1,alpha=0.025,beta=0.1){
  ## see friede and schmiedli 2010 (StatMed)
  nc <- (1+k*theta)^2 * (qnorm(alpha,lower=F)+qnorm(beta,lower=F))^2 / (t*lambda*(k+1)*theta*log(theta)^2)
  return(ceiling(c(nc=nc,nt=k*nc)))
}

condpowerrule_counts <- function(x,y,theta=1.5,maxN=Inf){
  k <- length(y)/length(x)
  lambda  <- mean(y) + mean(x)
  samplesize_counts(theta,lambda,k)
}
