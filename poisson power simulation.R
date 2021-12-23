
rm(list=ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is the statement:
# At least 44 subject is required for the study to have approximately 90% power to detect at least 25% reduction 
# in rates of events. # With a mean baseline event rate of 10 and sd of difference of 5 by 26 weeks 
# [ what is this sd of differences quoted for counts?] (25% is conservative, reduction expected to be bigger)
# 
# calc based on Wilcoxon signed rank test with 5K simulation (from Poisson) 2 sided alpha=5%
# robustness was also validated, using a one sample t test  
# 
# % reduction from baseline will be tested Hodges Lehman estimate of median percent reduction will also be reported

# note discrepancy above stating 2 sided and hypothesis below stated as 1 sided
# Ho p <.25
# H1 p >=.25
# Wilcoxon p value and >=25% is success
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
  
  # create a data frame
  my_data <- data.frame( 
    group = rep(c("before", "after"), each = n),
    counts = c(pp1,  pp2),
    ID=rep(1:n,2)
  )
  
  # analyse
  # https://stats.stackexchange.com/questions/71194/fitting-a-poisson-distribution-with-lme4-and-nlme
  # https://stats.stackexchange.com/questions/27869/fitting-a-poisson-glm-mixed-model-with-a-random-slope-and-intercept
  
  library(lme4)
  A <- glmer(counts ~ group + (1|ID), data=my_data, family="poisson")
  B <- glmer(counts ~ 1     + (1|ID), data=my_data, family="poisson")
  p <- anova(A,B)
  mix <-  p$`Pr(>Chisq)`[2]
  
  x <- summary(A)
  rate.reduction <- 1-1/exp(x$coefficients[2,1])
  
  x1 <- pp2 - pp1                                     # post - pre
  
  wil <- wilcox.test(x=x1, paired=FALSE)$p.value
  
  t.t <- t.test(x=x1, paired=FALSE)$p.value 
  
  newList <- list("glmer" = mix , "Wilcoxon signed rank test" = wil, "t.test"=t.t, "rate reduction"= rate.reduction)
  return(newList)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# execute
 
  res <- replicate(sims, 
                   po.power(n=22, r=.5, mu0=10, mu1=.75 ) )  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read outs
  alpha <- 0.05
  x <- NULL
  x <- t((res))
  x <- as.data.frame(x)
  mean( unlist(x[,"rate reduction"]))                      # mean rate reduction (%)
  median( unlist(x[,"Wilcoxon signed rank test"]))         # median Wilcoxon p-value
  table( unlist(x[,"glmer"])< alpha)/sims                     # power
  table( unlist(x[,"Wilcoxon signed rank test"])<alpha)/sims # power
  table( unlist(x[,"t.test"])<alpha)/sims                    # power

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 ## repeat but using simstudy package to generate correlated poisson

  po.power2 <- function(n=22, r=.75, mu0=10, mu1=.75) { # r correlation mu0 placebo rate, mu1 expected change in rate
    
    L1 <- mu0
    L2 <- mu0*mu1
    
    # l <- c(L1,L2)
    # 
    # my_data <- genCorGen(n*2, nvars = 2, params1 = l, dist = "poisson", rho = r, corstr = "cs", wide = FALSE)
    # names(my_data) <- c("ID","group", "counts")
    
    
    GenerateMultivariatePoisson<-function(p, samples, R, lambda){
      normal_mu=rep(0, p)
      normal = mvrnorm(samples, normal_mu, R)
      unif=pnorm(normal)
      pois=t(qpois(t(unif), lambda))
      return(pois)
    }
    
    
    d <-  GenerateMultivariatePoisson(p=2, samples=n, R=matrix(c(1, r, r, 1), 2, 2), lambda=c(L1,L2))
    
    my_data <- data.frame( 
      group = rep(c("before", "after"), each = n),
      counts = c(d[,1],  d[,2]),
      ID=rep(1:n,2)
    )
    
    
    library(lme4)
    A <- glmer(counts ~ group + (1|ID), data=my_data, family="poisson")
    B <- glmer(counts ~ 1     + (1|ID), data=my_data, family="poisson")
    p <- anova(A,B)
    mix <-  p$`Pr(>Chisq)`[2]
    
    x <- summary(A)
    rate.reduction <- 1-1/exp(x$coefficients[2,1])
    
    x1 <- d[,1] - d[,2]                                     # post - pre
    
    wil <- wilcox.test(x=x1, paired=FALSE)$p.value
    
    t.t <- t.test(x=x1, paired=FALSE)$p.value 
    
    newList <- list("glmer" = mix , "Wilcoxon signed rank test" = wil, "t.test"=t.t, "rate reduction"= rate.reduction)
    return(newList)
    
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # execute
   
  res <- replicate(sims, 
                   po.power2(n=22, r=.5, mu0=10, mu1=.75 ) )  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # read outs
  alpha <- 0.05
  x <- NULL
  x <- t((res))
  x <- as.data.frame(x)
  mean( unlist(x[,"rate reduction"]))                      # mean rate reduction (%)
  median( unlist(x[,"Wilcoxon signed rank test"]))         # median Wilcoxon p-value
  table( unlist(x[,"glmer"])< alpha)/sims                     # power
  table( unlist(x[,"Wilcoxon signed rank test"])<alpha)/sims # power
  table( unlist(x[,"t.test"])<alpha)/sims                    # power
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

























#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this is not a REAL paired simulation!! the responses are independent - so dont use
#      quick.po <- function() {
#       y_0 <- rpois(n = 22, lambda = 10)
#       y_1 <- rpois(n = 22, lambda = 7.5)
#       
#     A <- wilcox.test(y_0, y_1, paired=TRUE)$p.value 
#     B <- t.test(     y_0, y_1, paired=TRUE)$p.value
#       
#       #A <- wilcox.test(y_0-y_1, paired=FALSE)$p.value 
#       #B <- t.test(     y_0- y_1, paired=FALSE)$p.value
#       
#       #A <- wilcox.test(y_0/y_1, paired=FALSE)$p.value 
#       #B <- t.test(     y_0/y_1, paired=FALSE)$p.value
#     
#       newList <- list(  "Wilcoxon signed rank test" = A, "t.test"=B)
#       return(newList)
#    }
#      
#      sims <- 999
#      res <- replicate(sims, 
#                       quick.po  ( ) )  
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      # power read outs
#      x <- t(data.frame(res))
#      x <- (x<0.05)
#      x <- data.frame(x)
#      cor(x)
#      
#       
#      table(x$ Wilcoxon.signed.rank.test)/sims
#      table(x$ t.test)/sims
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      
#      
#       
#       
#       
#       simpois  <- function(obs=44, mu=10, rate=.75)
#    
#       lambda_baseline <- mu 
#       y_0 <- rpois(n = obs, lambda = lambda_baseline)
#       ## Simulate the group allocation
#       group <- rep(1,obs) #c(rep(1, obs/2), rep(0, obs/2))
#       beta <- log(rate)
#       ## Simulate the follow-up counts
#       lambda_followup <- exp(beta * group) * mu 
#       y_1 <- rpois(n = obs, lambda = lambda_followup) 
#   
#       summary( glm(y_1 ~ y_0, family="poisson"))
      
      
      
  
      
      
      
      # obs <- 44
      # rate1 <- 10
      # y_0 <- rpois(n = obs, lambda = rate1)
      # ## Simulate the group allocation
      # group <- 1 # c(rep(1, obs/2), rep(0, obs/2))
      # beta <- log(7.5)
      # ## Simulate the follow-up counts
      # lambda_followup <- exp(beta * group)  
      # y_1 <- rpois(n = obs, lambda = lambda_followup) 
      # 
      # mean(y_0)
      # mean(y_1)
      # 
      # summary(MASS::glm.nb(y_1 ~ y_0))
      # 
      # r <- y_0 / y_1 
      # summary(glm.nb( r ~ 1))
      # 
      # summary( glm(r ~ 1, family="poisson"))
      # summary( glm(y_1 ~ y_0, family="poisson"))
      # 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # is this the calculation
      
      # n     <- 22
      # sdd   <- 5
      # delta <- 10-7.5
      # 
      # mean(replicate(10000, t.test(     x=rnorm(n, delta, sd), paired=FALSE)$p.value <0.05    ))
      # mean(replicate(10000, wilcox.test(x=rnorm(n, delta, sd), paired=FALSE)$p.value <0.05    ))
      # 
      # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # # this was another study I qcd erik
      # 
      # n     <- 154
      # delta <- 8.5
      # sd    <- 3.14
      # sdd   <- sd*sqrt(2)
      # 
      # mean(replicate(10000, diff(t.test(x=rnorm(n, delta, sd),  paired=FALSE)$conf.int)))
      # mean(replicate(10000, diff(t.test(x=rnorm(n, delta, sdd), paired=FALSE)$conf.int)))
      # 
      # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 
      # # simulate correlated counts library(simstudy)
      # 
      # 
      # l <- c(10, 7.5) #
      # dx <- genCorGen(44, nvars = 2, params1 = l, dist = "poisson", rho = .3, corstr = "cs", wide = TRUE)
      # dx
      # dx <- genCorGen(44, nvars = 2, params1 = l, dist = "poisson", rho = .3, corstr = "cs", wide = FALSE)
      # dx
      # 
      
      