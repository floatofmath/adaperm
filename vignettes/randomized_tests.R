library(adaperm)
g1 <- rep(0:1,each=5)
g2 <- rep(0:1,each=5)
gE <- rep(0:1,each=3)
x1 <- rnorm(10)
x2 <- rnorm(10)
xE <- rnorm(6)

t.test(-c(x1,x2,xE)~c(g1,g2,gE))
permutations <- 10^6
restricted=T
stat <- meandiff

pdist  <- perm_dist(x1,x2,g1,g2,stat,B=permutations,restricted=restricted)
## conditional distribution of preplanned test given the first stage data
cdist <- adaperm:::cond_dist(x1,x2,g1,g2,stat,B=permutations,restricted=restricted)
## permutation distribution of adapted test
edist <- perm_dist(x2,xE,g2,gE,stat,B=permutations,restricted=restricted)
## observed test statistic of adapted test
t2 <- stat(c(x2,xE),c(g2,gE))

## adaptive p-value non-randomized test
pval <- mean(pdist >= quantile(cdist,mean(edist<t2),type=1))

p2 <- 1-mean(edist<t2)
quantile(cdist,mean(edist<t2),type=1)
mean(cdist>= quantile(cdist,mean(edist<t2),type=1))
ps = mean(pdist >= max(cdist[rank(-cdist,ties='max') > mean(edist>=t2)*length(cdist)]))

Tas <- max(cdist[rank(-cdist,ties='max') > mean(edist>=t2)*length(cdist)])
## possible cer values
cers <- rank(-cdist,ties='max')/length(cdist)


amax <- max(cers[cers < mean(edist>=t2)])
Tmax <- min(cdist[cers < mean(edist>=t2)])
as <- min(cers[cers >= mean(edist>=t2)])
Tas <- max(cdist[cers >= mean(edist>=t2)])
pmax <- mean(pdist >  Tmax) + trim
trim <- (mean(edist>=t2) - amax)*mean(cdist==Tmax)/mean(cdist==Tmax)
## stepping down one step on the CER scale, means stepping down quite af few on the p scale,
## the difference can not be overcome by randomization
inbetween <- sort(pdist[pdist <= Tmax & pdist >= Tas])
## if we where to perform a test at a couple of steps lower level, we are not able to reject
mean(cdist >= Tas) > p2
mean(cdist >= Tmax) > p2
sapply(inbetween,function(t) mean(cdist >= t) > p2)
inbetween == Tas

#########################################################################
### Discrete data
#########################################################################
g1 <- rep(0:1,each=5)
g2 <- rep(0:1,each=5)
gE <- rep(0:1,each=3)
x1 <- sample(-1:1,10,rep=T)
x2 <- sample(-1:1,10,rep=T)
xE <- sample(-1:1,6,rep=T)

t.test(-c(x1,x2,xE)~c(g1,g2,gE))
permutations <- 10^6
restricted=T
stat <- meandiff

pdist  <- perm_dist(x1,x2,g1,g2,stat,B=permutations,restricted=restricted)
## conditional distribution of preplanned test given the first stage data
cdist <- adaperm:::cond_dist(x1,x2,g1,g2,stat,B=permutations,restricted=restricted)
## permutation distribution of adapted test
edist <- perm_dist(x2,xE,g2,gE,stat,B=permutations,restricted=restricted)
## observed test statistic of adapted test
t2 <- stat(c(x2,xE),c(g2,gE))

## adaptive p-value non-randomized test
pval <- mean(pdist >= quantile(cdist,mean(edist<t2),type=1))

p2 <- 1-mean(edist<t2)
quantile(cdist,mean(edist<t2),type=1)
mean(cdist>= quantile(cdist,mean(edist<t2),type=1))
ps = mean(pdist >= max(cdist[rank(-cdist,ties='max') > mean(edist>=t2)*length(cdist)]))

Tas <- max(cdist[rank(-cdist,ties='max') > mean(edist>=t2)*length(cdist)])
## possible cer values
cers <- rank(-cdist,ties='max')/length(cdist)


amax <- max(cers[cers < mean(edist>=t2)])
Tmax <- min(cdist[cers < mean(edist>=t2)])
as <- min(cers[cers >= mean(edist>=t2)])
Tas <- max(cdist[cers >= mean(edist>=t2)])
trim <- (mean(edist>=t2) - amax)*mean(cdist==Tmax)/mean(cdist==Tmax)
pmax <- mean(pdist >  Tmax) + trim

## stepping down one step on the CER scale, means stepping down quite af few on the p scale,
## the difference can not be overcome by randomization
inbetween <- unique(sort(pdist[pdist <= Tmax & pdist >= Tas]))
## if we where to perform a test at a couple of steps lower level, we are not able to reject
mean(cdist >= Tas) > p2
mean(cdist >= Tmax) > p2
sapply(inbetween,function(t) mean(cdist >= t) > p2)
inbetween == Tas

## one possibility would be to look for the cer that yields the mid-p value

-4       -2    0    2 4
15/16 11/16 5/16 1/16 0

a <- c(-4,-2,-2,-2,-2,0,0, 0, 0, 0, 0, 2, 2, 2, 2, 4)
b <- c(1, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6,12,12,12,12,16)
c <- c(1, 5, 5, 5, 5,11, 11,11,11,11,11,15,15,15,15,16)

rank(-a,ties='max')
rank(-a,ties='min')

Ta <- function(p) max(a[rank(-a,ties='max') > p*length(a)])
Tas <- function(p) min(a[rank(-a,ties='min') < p*length(a)])

rank(-a,ties='max')/length(a)
p=1/16
Ta(1/17)
Tas(1/17)

Ta(1/16)
Ta(1/16)
Ta(2/16)
Tas(2/16)
Ta(5/16)
Ta(5/16)
Ta(6/16)
Tas(6/16)


