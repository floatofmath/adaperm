# ##### alternative version of the function!!
# adaptive_permtest_2s <- function(x,y,n1,n,ne,test_statistic,
#                                  m1=n1,m=n,me=ne,perms=1000,alpha=0.025,
#                                  restricted, 
#                                  type ='non-randomized'){
#   if(ne>n){
#     xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
#   } else {
#     if(me>m) stop('Stages with controls only not supported')
#     xs <- split(x,rep(1:2,c(n1,ne-n1)))
#   }
#   if(me>m){
#     ys <- split(y,rep(1:3,c(m1,m-m1,me-m)))
#   } else {
#     if(ne>n) stop('Stages with treatments only not supported')
#     ys <- split(y,rep(1:2,c(m1,me-m1)))
#   }
#   gs <- lapply(1:length(xs),function(i) rep(0:1,c(length(xs[[i]]),length(ys[[i]]))))
#   xs <- lapply(1:length(xs),function(i) c(xs[[i]],ys[[i]]))
#   if(me==m){
#     return(alpha>=perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,
#                             B=perms,restricted=restricted,type=type))
#   }
#   A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,permutations=perms,alpha=alpha)
#   q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,B=perms,restricted=restricted,type=type)
#   A>=q
# }
# ##########################
# 
# adaptive_permtest2_2s <- function(x,y,n1,n,ne,test_statistic,
#                                   m1=n1,m=n,me=ne,perms=1000,alpha=0.025,
#                                   resam=100,restricted=FALSE,
#                                   type=c('non-randomized'),
#                                   atest_type='midp'){
#   if(ne>n){
#     xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
#   } else {
#     if(me>m) stop('Stages with controls only not supported')
#     xs <- split(x,rep(1:2,c(n1,ne-n1)))
#   }
#   if(me>m){
#     ys <- split(y,rep(1:3,c(m1,m-m1,me-m)))
#   } else {
#     if(ne>n) stop('Stages with treatments only not supported')
#     ys <- split(y,rep(1:2,c(m1,me-m1)))
#   }
#   gs <- lapply(1:length(xs),function(i) rep(0:1,c(length(xs[[i]]),length(ys[[i]]))))
#   xs <- lapply(1:length(xs),function(i) c(xs[[i]],ys[[i]]))
#   if(me==m){
#     return(alpha>=perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,B=perms,restricted=restricted,type=type))
#   }
#   A <- permutation_CER2(xs[[1]],gs[[1]],xs[[2]],xs[[3]],test_statistic,one_sample=FALSE,restricted=FALSE,
#                         permutations=perms,subsamples=resam,alpha=alpha)
#   q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,B=perms,restricted=restricted,type=type)
#   A>=q
# }


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
                                      m1=n1,m=n,resam=100, perms,...){
  x <- do.call(match.fun(rdist),c(n=n,control_opts))
  y <- do.call(match.fun(rdist),c(n=n,treatment_opts))
  nes <- rule(x[1:n1],y[1:m1])
  print(nes)
  if(nes[1]>n){
    ne <- nes[1]
    x <- c(x,do.call(match.fun(rdist),c(n=nes[1]-n,control_opts)))
  } else {
    ne <- n
  }
  if(nes[2]>m){
    me <- nes[2]
    y <- c(y,do.call(match.fun(rdist),c(n=nes[2]-n,treatment_opts)))
  } else {
    me <- m
  }

  n_combs <- prod(choose(c(n1+m1,n-n1+m-m1,ne+me),c(n1,n-n1,ne)))
  allperms <- (n_combs > perms) #### resam??
  if(allperms){
    Aother <- 'midp'
    cer_type <- "randomized"
      } else {
      Aother <- 'davison_hinkley'
      cer_type <- "non-randomized"
      }

  list(ne = ne,
       perm_CERnonrndmz_Anonrndmz = adaperm_DR(c(x,y),c(rep(0,me),rep(1,ne)),n1=n1,n=n,m1=m1,m=m,
                                            test=meandiff,cer_type='non-randomized',atest_type="non-randomized",
                                            permutations=perms),
       perm_CERother_Aother = adaperm_DR(c(x,y),c(rep(0,me),rep(1,ne)),n1=n1,n=n,m1=m1,m=m,
                                            test=meandiff,cer_type='randomized',atest_type=Aother,
                                            permutations=perms),
       perm_ratio_CERnonrndmz_Anonrndmz =  adaperm_DR(c(x,y),c(rep(0,me),rep(1,ne)),n1=n1,n=n,m1=m1,m=m,
                                                   test=meanratio,cer_type='non-randomized',atest_type="non-randomized",
                                                   permutations=perms),
       perm_ratio_CERother_Aother =  adaperm_DR(c(x,y),c(rep(0,me),rep(1,ne)),n1=n1,n=n,m1=m1,m=m,
                                                   test=meanratio,cer_type='randomized',atest_type=Aother,
                                                   permutations=perms),
       invnorm = adaptive_invnormtest_2s(x,y,n1,n,ne,m1,m,me),
       invnorm_wlcx = adaptive_invnorm_wilcoxtest_2s(x,y,n1,n,ne,m1,m,me),
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
  ## see Schneider + schmiedli + friede 2013 (biometrical journal)
  ## n 0 = (z 1−β + z 1−α ) ^2 /   log(θ ∗ ) ^2 *(1 + kθ ∗)  /θ ∗ k/tλ ∗ 0
  nc <-(qnorm(alpha,lower=F)+qnorm(beta,lower=F))^2 *(1+k*theta) / (t*lambda*k*theta*log(theta)^2)
  return(ceiling(c(nc=nc,nt=k*nc)))
}

condpowerrule_counts <- function(x,y,theta=1.5,maxN=Inf){
  k <- length(y)/length(x)
  lambda  <- mean(x)
  pmin(maxN,samplesize_counts(theta,lambda,k))
}


