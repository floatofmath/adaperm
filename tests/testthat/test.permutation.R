context("Permutation Spaces")
## test omega
n1 <- 9
g1 <- sample(c(1,-1),n1,rep=T)

test_that("All (re)assignments",
          {
              k <- 4
              o <- all_reassignments(n1,k)
              uo <- all_assignments(n1)
              expect_equal(dim(o),c(n1,choose(n1,k)))
              expect_equal(k,unique(colSums(o)))
              expect_equal(dim(uo),c(n1,2^n1))
              expect_less_than(max(colSums(uo)),n1+1)
              B <- 100
              ro <- random_reassignments(n1,k,B)
              expect_equal(dim(ro),c(n1,B))
              expect_equal(k,unique(colSums(ro)))
              ruo <- random_assignments(n1,B)
              expect_equal(dim(ruo),c(n1,B))
          })


test_that("Stratified (re)assignment",
          {
              ns <- c(5,4,3,2)
              ks <- c(3,2,1,1)
              o <- strat_reassignments(ns,ks)
              expect_equal(dim(o),c(sum(ns),prod(choose(ns,ks))))
              expect_equal(sum(ks),unique(colSums(o)))
              expect_equal(ks,sapply(list(colSums(o[1:5,]),colSums(o[6:9,]),colSums(o[10:12,]),colSums(o[13:14,])),unique))
          })

test_that("Omega one stage, restricted",
          {
              o1 <- omega(g1)
              expect_equal(sum(g1>0),unique(colSums(o1 > 0)))
              expect_equal(ncol(o1),choose(n1,sum(g1>0)))
              expect_equal(nrow(o1),length(g1))
          })


test_that("Omega one stage unrestricted",
          {
              B <- 11
              o2 <- omega(g1,restricted=FALSE)
              expect_equal(ncol(o2),2^length(g1))
              expect_equal(nrow(o2),length(g1))
              o2 <- omega(g1,restricted=FALSE,B=10^6)
              expect_equal(dim(o2),c(length(g1),2^length(g1)))
              o2 <- omega(g1,restricted=FALSE,B=B)
              expect_equal(dim(o2),c(length(g1),B+1))
          })


n2 <- 7
g2 <- sample(c(1,-1),n2,rep=T)

test_that("Omega two stages restricted",
          {
              o3 <- omega(g1,g2,B = 10^6)
              expect_equal(sum(g1>0) + sum(g2>0),unique(colSums(o3 >0)))
              expect_equal(dim(o3),c(length(g2) + length(g1),choose(n1,sum(g1>0))*choose(n2,sum(g2>0))))
          })
## test conddist:
pasta  <- function(y,xs) apply(xs,2,function(x) paste(which(x>0),collapse=''))
test_that("Conditional permutation distribution",
          {
              x1 <- rnorm(4)
              x2 <- rnorm(4)
              expect_true(all(cond_dist(x1,x2,rep(c(-1,1),each=2),rep(c(-1,1),each=2),pasta,100000) %in% apply(expand.grid('34',apply(4+gtools::combinations(4,2),1,paste,collapse='')),1,paste,collapse='')))
          })

test_that("Permutation distribution",
          {
              expect_true(all(perm_dist(rnorm(4),rnorm(4),rep(c(-1,1),each=2),rep(c(-1,1),each=2),pasta,100000) %in% apply(expand.grid(apply(gtools::combinations(4,2),1,paste,collapse=''),apply(4+gtools::combinations(4,2),1,paste,collapse='')),1,paste,collapse='')))
          }
          )

n <- 8
n1 <- 4
x1 <- rep(c(1,-1),n1/2)
g1 <- rep(c(-1,1),each=n1/2)
x2 <- rep(c(1,-1),n1/2)
g2 <- rep(c(-1,1),each=n1/2)
test_that("Permutation test",
          {
              expect_equivalent(perm_test(x1,x2,g1,g2,sumdiff,100000,restricted=T),.75)
              expect_equivalent(perm_test(x1,x2,g1,g2,sumdiff,100000,restricted=T,type='midp'),.5)
              expect_equivalent(perm_test(x1,x2,g1,g2,sumdiff,100000,restricted=T,type='davison_hinkley'),10/19)
          })
          
## test t_test
n <- 8
n1 <- 4
x1 <- rep(c(1,-1),n1/2)
g1 <- rep(c(-1,1),each=n1/2)
x2 <- rep(c(1,-1),n1/2)
g2 <- rep(c(-1,1),each=n1/2)
test_that("t-test",
          {
              expect_equivalent(t_test(x1,x2,g1,g2),.5)
          }
          )



## test permutation cer
n <- 12
n1 <- 6
x1 <- rep(c(-1,1),each=n1/2)
g1 <- rep(c(-1,1),each=n1/2)
x2 <- rep(c(-1,1),each=n1/2)
g2 <- rep(c(-1,1),each=n1/2)
test_that("permutation CER",
          {
              expect_equivalent(permutation_cer(x1,x2,g1,test_statistic=sumdiff,permutations=100000,alpha=.025,restricted=T),0.05)
              expect_equivalent(permutation_cer(x1,x2,g1,test_statistic=sumdiff,permutations=100000,alpha=.025,restricted=T,cer_type='randomized'),0.275)
              expect_equivalent(adaperm_DR(c(x1,x2),g1,n1=3,n=6,test_statistic=sumdiff,atest_type='CER'),0.05)
          })


## test t2p
pd <- perm_dist(x1,x2,g1,g2,sumdiff,1000)
qs <- runif(1000)
ta <- sapply(qs,p2t,dist=pd)
pt <- sapply(ta+0.001,t2p,dist=pd,type='davison_hinkley')
test_that("permutation p-value and T_alpha",
          {
              expect_true(all(pt<=qs))
          })
          
