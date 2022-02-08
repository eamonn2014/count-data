#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  generate cor data , see references
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# given R, calc sd of differences
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# R is the correlation of X Y

#var(x - y) = var(x) + var(y) - 2 x [R x sd(x) x sd(y)]  [1] # in [] equals covariance

# Q. get correlation in 1 ?

#corr(x,y) = cov(x,y) / sd(x) * sd(y)  [2]

# rearrange [2]
#cov(x,y) = corr(x,y) * sd(x) * sd(y)  [2a]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# so substitute [2a] into end of 1

#var(x - y) = var(x) + var(y) - 2 x R x sd(x) x sd(y)  [3]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# so sd of diff is sqrt of [3]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# refs:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#[1] https://stats.libretexts.org/Bookshelves/Introductory_Statistics/Book%3A_Introductory_Statistics_(Lane)/04%3A_Describing_Bivariate_Data/4.07%3A_Variance_Sum_Law_II_-_Correlated_Variables

#[2] https://www.countbayesie.com/blog/2015/2/21/variance-co-variance-and-correlation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(lme4)

  sims.p=100
  n     =22
  r     =0.75
  mu0.p =10
  eff.p =0.75 
  alpha =0.05


formatz3 <- function(x){
  sprintf(x, fmt = '%#.3f')  
}



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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to generate correlated poisson and analyse and evaluate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

po.cor.power <- function(n=22, r=.75, mu0=10, eff.p=.75) { # r correlation mu0 placebo rate, mu1 expected change in rate
  
  #n=22; r=.75; mu0=10; eff.p=.75
  # 1. create correlated poisson
  # 2. analyse many times and examine
  
  L1 <- mu0
  L2 <- mu0*eff.p
  
  # generate correlated poisson ref: https://thomasward.com/simulating-correlated-data/
  # Sample correlated N(0, 1) distributions from a multivariate normal distribution.
  # Transform them to correlated Uniform(0, 1) distributions with the normal CDF.
  # Transform them to any correlated probability distribution you desire with that probability distribution’s inverse CDF.
  
  Sigma <- matrix(c(1, r, r, 1), 2, 2)
  
  p2 <- mvrnorm(n, Sigma = Sigma)   # correlated continuous, see function in global area!
  
  U <- pnorm(p2, mean = 0, sd = 1)  # correlated uniform
  pp1 <- qpois(U[, 1], L1)          # correlated poisson
  pp2 <- qpois(U[, 2], L2)          # correlated poisson
  
  cor. <- cor(pp1,pp2)              # capture correlation
  
  # create a data frame
  my_data <- data.frame( 
    group = rep(c("A.before", "B.after"), each = n),
    counts = c(pp1,  pp2),
    ID=rep(1:n,2)
  )
  
  # mixed model
  A <- glmer(counts ~ group + (1|ID), data=my_data, family="poisson")
  
  # B <-    glmer(counts ~ 1     + (1|ID), data=my_data, family="poisson")  # dont use LRT to speed up sims
  # p <-    anova(A,B)
  # mix <-  p$`Pr(>Chisq)`[2]
  
  x <- summary(A)
  mix <- summary(A)$coeff[2,"Pr(>|z|)"]  # to speed up simulations we don't use LR test but Wald test
  
  rate.reduction <- exp(x$coefficients[2,1])
  
  intercept <- exp( x$coefficients[1,1])
  beta      <- exp( x$coefficients[2,1])
  
  x1 <- pp2 - pp1                                                                         # post - pre
  
  w <- wilcox.test(x=x1, paired=FALSE, correct=FALSE, conf.int=TRUE  ,conf.level = 0.95)  # paired is false as we have taken the difference in pairs
  
  wil <- w$p.value
  wile <- w$estimate[1][[1]]
  
  t.t <- t.test(x=x1, paired=FALSE)$p.value # t.test(x1)$p.value 
  
  signed_rank = function(x) sign(x) * rank(abs(x))
  ttor <- t.test(signed_rank(x1))$p.value
  
  # wilcox.test(x1)
  # t.test(signed_rank(x1))
  
  #lets capture SD of differences
  sdd <- sd(x1)
  
  #newList <- list("glmer" = mix , "Wilcoxon signed rank test" = wil,
  #    "t.test"=t.t, "rate reduction"= rate.reduction, "intercept"=intercept, "beta"=beta,
  #   "ttest on ranks"=ttor, "sd of diff"=sdd)
  #return(newList)
  # collect all relevant info
  c(rate.reduction,mix , wil,t.t,  ttor,intercept, beta ,sdd, wile, cor.)
  
}


  res <- replicate(sims.p, 
                   po.cor.power(n=n, r=r, mu0=mu0.p,  eff.p=eff.p ) )  
  
  x <- NULL
  x <- t((res))
  x <- as.data.frame(x)

  
  R <- r             # seq(-.9,.9,.1)   # correlation
  X <- mu0.p         # var of poisson (= mean)
  Y <- mu0.p*eff.p   # var of poisson (= mean)
  
  cat(paste0("Evaluations based on ",sims.p," simulations take into consideration correlation, the default starting value is "
             ,r,". The observed mean correlation is ", formatz3(mean( unlist(x[,"V10"]))), "\nThe median p-value from a mixed effects Poisson model is "
             ,(median( unlist(x[,"V2"])))   , 
             ". \nPower based on the mixed model Wald test is ",
             
             formatz3((table( unlist(x[,"V2"])< alpha)/sims.p)[2][[1]]),
             
             "\nPower based on the Wilcoxon signed rank test power is ",
             formatz3((table( unlist(x[,"V3"])< alpha)/sims.p)[2][[1]]),
             "\nMean Hodges-Lehmann estimator (pseudo-median) [1) Calc. post-pre; 2) Calc. Walsh averages, i.e. n*(n+1)/2 pairwise averages; 3) report median] ",
             formatz3(mean( unlist(x[,"V9"])))," with 95%CI ",formatz3(quantile(unlist(x[,"V9"]), .025))," to ",formatz3(quantile(unlist(x[,"V9"]), .975)),      "\nPower based on the t-test is ",
             formatz3((table( unlist(x[,"V4"])< alpha)/sims.p)[2][[1]]),
             "\nPower based on the t-test of ranks is ",
             formatz3((table( unlist(x[,"V5"])< alpha)/sims.p)[2][[1]]),
             "\nThe mean intercept from the mixed effects Poisson model is ",
             formatz3(mean( unlist(x[,"V6"]))),
             "\nThe mean treatement effect from the mixed effects Poisson model is ",
             formatz3(mean( unlist(x[,"V7"])))," with 95%CI ",formatz3(quantile(unlist(x[,"V7"]), .025))," to ",formatz3(quantile(unlist(x[,"V7"]), .975)),
             "\nMean SD of differences is ", 
             formatz3(mean( unlist(x[,"V8"]))),  
             "\nThe theoretical SD of differences based on Poisson variances and stated correlation is ",
             formatz3(sqrt(X+Y-(2*sqrt(X)*sqrt(Y)*R))),
             "\n"))
 