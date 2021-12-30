
# here we are looking at the tandard deviation of the difference
# Var(X−Y)=Var(X)+Var(Y)−2Cov(X,Y)
# and show that we ~match the theory

R <- seq(-.9,.9,.1)   # correlation
X <- 10             # var of poisson (=mean)
Y <- 7.5            # var of poisson (=mean)

# correlation = covariance/ sd(x)*sd(y)*covariance, so covar = cor*sd(x)*sd(y)
# expected sd of differences

sqrt(X+Y-(2*sqrt(X)*sqrt(Y)*R))



rm(list=ls())
library(MASS)
library(lme4)
sims=999
po.power <- function(n=22, r=.75, mu0=10, mu1=.75) { # r correlation mu0 placebo rate, mu1 expected change in rate
  
  # 1. create correlated poisson
  # 2. analyse many times and examine
  
  L1 <- mu0
  L2 <- mu0*mu1
  
  # generate correlated poisson ref: https://thomasward.com/simulating-correlated-data/
  
  # Sample correlated N(0, 1) distributions from a multivariate normal distribution.
  # Transform them to correlated Uniform(0, 1) distributions with the normal CDF.
  # Transform them to any correlated probability distribution you desire with that probability distribution’s inverse CDF.
  
  Sigma <- matrix(c(1, r, r, 1), 2, 2)
  #Sigma
  
  mvrnorm <- function(n = 1, mu = 0, Sigma) {
    nvars <- nrow(Sigma)
    # nvars x n matrix of Normal(0, 1)
    nmls <- matrix(rnorm(n * nvars), nrow = nvars)
    # scale and correlate Normal(0, 1), "nmls", to Normal(0, Sigma) by matrix mult
    # with lower triangular of cholesky decomp of covariance matrix
    scaled_correlated_nmls <- t(chol(Sigma)) %*% nmls
    # shift to center around mus to get goal: Normal(mu, Sigma)
    samples <- mu + scaled_correlated_nmls
    # transpose so each variable is a column, not
    # a row, to match what MASS::mvrnorm() returns
    t(samples)
  }
  
  p2 <- mvrnorm(n, Sigma = Sigma)   # correlated continuous
  
  U <- pnorm(p2, mean = 0, sd = 1)  # correlated uniform
  
  pp1 <- qpois(U[, 1], L1)          # correlated poisson
  pp2 <- qpois(U[, 2], L2)          # correlated poisson
   
  x1 <- pp2 - pp1                                     # post - pre
  
  mu= mean(x1)
  sdd <- sd(x1)
  cor. <- cor(pp1, pp2)
  
#  newList <- list( "mean diff"=mu, "sd of diff"=sdd, "cor"=cor.)
 # return(newList)
  c(mu,sdd)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# execute

res <- replicate(sims, 
                 po.power(n=220, r=-.4, mu0=10, mu1=.75 ) )  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read outs
alpha <- 0.05
x <- NULL
x <- t((res))
x <- as.data.frame(x)
 
summary( unlist(x[,"mean diff"])) 
summary( unlist(x[,"sd of diff"])) 
summary( unlist(x[,"cor"])) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# do a bigger simulation

R <- seq(-.9,.9,.1)
reps <- 1:3
mu0 <-.75*mu1


dnam = list( R=R, reps=reps , What=c("sd of diff","mean diff") )
pwpr <- array( NA, dim=sapply( dnam, length ), dimnames=dnam )
str( pwpr )



system.time(
  for( i in 1:length(R) )
    for( j in 1:length(reps) )
        pwpr[i,j,] <- po.power( n=220, r=R[i],  mu0=10, mu1=.75 )
                  )
              
L <- reshape::melt(pwpr)
with(L, tapply(value,list(R, What), mean))
 














