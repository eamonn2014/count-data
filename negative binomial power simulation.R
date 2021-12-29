#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# make a function for negative binomial power
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Patient heterogeneity beyond that captured by patient level covariates can be expected

# Simulate a RCT trial, randomisation 1:1, for a count endpoint. We will use the 
# negative binomial model with natural log of duration of follow up time as an offset variable
# the duration of follow up will be the planned duration of follow up for all patients.
# The offset will make regression coefficients easier to interpret.

# n number of patients per arm (needs to be even for this simulation)
# dispersion parameter, believe we should talk about 1 over this value?
# mu0 rate for placebo - events per patient year of follow up
# mu1 rate reduction compared to placebo that is hypothesised for the new treatment.
# drop out rates for treatment discontinuation
# follow up period (in years)
# what do we estimate for mu1? This is the relative risk ratio!
library(MASS)
nb.power <- function(n=220, disp=1.3, mu0=1, mu1=.65, drop1=.1, drop2=.1, fup=1) {
  
 # n=220; disp=1.3; mu0=1; mu1=.65; drop1=.001; drop2=.001; fup=1
  
  dose <- c(rep("placebo",n),rep("trt",n)) # 50:50 split of patients
  
  mu   <- c(rep(mu0,n), rep(mu0*mu1,n))    # rates in two arms
  
  drop <- c(rep(drop1,n), rep(drop2,n))    # tr discontinuation rates
  
  # constant exponential hazard rate over time for treatment discontinuation is assumed with drop1 and 
  # drop2 %s of patients off treatment by 1 year in active and placebo. 
  # No further treatment effect versus placebo is assumed during off treatment, 
  # that is mu0 events per patient year is assumed for patients after treatment discontinuation in both arms.
  
  f <- - rexp(2*n) / log(1-drop)
  length <- ifelse(f > fup, fup, f)   # curtail at follow up time 
  
  # create responses for each patient
  # accounting for on treatment and time for those off treatment, notice the parameterization
  # size is the dispersion parameter
  y <-  rnbinom(n*2,  prob=1/(1+ mu*length* disp),        size=1/disp)  +   # on treatment period events
        rnbinom(n*2,  prob=1/(1+ mu0*(fup-length)*disp),  size=1/disp)      # accounting for any off treat. period events
  
  lfup <- log(fup)             # exposure log(fup) for everyone, that is fup year
  logleng  <- rep(lfup, n*2)   #  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~qc
  # L <- list(dose,mu,drop,length,logleng)
  # 
  # require(tidyverse)
  # par(mfrow=c(1,2))
  # y[1:n]     %>% table %>% barplot() #quick and dirty
  # y[n+1:n*2] %>% table %>% barplot() #quick and dirty
  # par(mfrow=c(1,1))
  # 
  # mean(y[1:n])  # matches the rate
  # mean(y[(n+1):(n*2)])
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # analyse with neg. binomial model
  x <- summary(MASS::glm.nb(y~dose+offset((logleng))))
  #exp(x$coeff)
  # collect p-values
  p <-  x$coefficients["dosetrt","Pr(>|z|)"]
  
  return(p)
  
}

res <- replicate(1000, nb.power(drop1=.1, drop2=.1))  #
mean(res<0.05)

# https://www.nejm.org/doi/pdf/10.1056/NEJMoa1403290?articleTools=true
# published power statement that we reproduce? The RCT in question is however 3 arm, but my simulation only looks at 2 arms

res <- replicate(999, nb.power(n=180, disp=1/.8, mu0=2.4, mu1=.6, drop1=0.0001, drop2=0.0001, fup=32/52) )  #
mean(res<0.05)

# We estimated that with 180 patients in each group, the study would 
# have a power of 90% to detect a 40% decrease in the exacerbation rate,
# from 2.40 per year in the placebo group to 1.44 per year in each of the mepolizumab groups, 
# at a two-sided significance level of 0.05. In performing this calculation, we assumed that the number of exacerbations 
# would follow a negative binomial distribution with a dispersion parameter k=0.8.
# the primary efficacy analysis was performed at 32 weeks as stated in the paper.
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end make a function
res <- replicate(1000, nb.power(mu0=1.5,n=120))  #.59
mean(res<0.05)
res <- replicate(1000, nb.power(mu0=.9,n=120))   #.49
mean(res<0.05)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3 arm simulation
# Dad's Army:
# Good heavens
# What do you mean good heavens? Good heavens


nb.power2 <- function(n=220, disp=1.3, mu0=1, mu1=.65, mu2=.65, drop1=.1, drop2=.1, drop3=.1, fup=1) {
  
  #n=220; disp=1.3; mu0=1; mu1=.65; mu2=.65; drop1=.1; drop2=.1; drop3=.1; fup=1
  #n=180; disp=1/.8; mu0=2.4; mu1=.6; mu2=.6; drop1=0.0001; drop2=0.0001; drop3=0.0001; fup=32/52
  
  dose <- c(rep("placebo",n), rep("trt1",n) , rep("trt2",n) )  # 1/3 split of patients
  
  mu   <- c( rep(mu0,n),  rep(mu0*mu1,n)  , rep(mu0*mu2,n))    # rates in two arms
  
  drop <- c( rep(drop1,n), rep(drop2,n), rep(drop3,n))         # trt discontinuation rates
  
  f <- - rexp(3*n) / log(1-drop)
  length <- ifelse(f > fup, fup, f)   # curtail at follow up time 
  
  
  y <-  rnbinom(n*3,  prob=1/(1+ mu*length* disp),        size=1/disp)  +   # on treatment period events
    rnbinom(n*3,  prob=1/(1+ mu0*(fup-length)*disp),  size=1/disp)      # accounting for any off treat. period events
  
  lfup <- log(fup)             # exposure log(fup) for everyone, that is fup year
  logleng  <- rep(lfup, n*3)   #  
  
  x <- summary(MASS::glm.nb(y~dose+offset((logleng))))
  
  # collect p-values
  p1 <-  x$coefficients["dosetrt1","Pr(>|z|)"]
  p2 <-  x$coefficients["dosetrt2","Pr(>|z|)"]
  
  newList <- list("p1" = p1 , "p2" = p2)
  return(newList)
  
}

sims <- 299
res <- replicate(sims, 
                 nb.power2(n=180, disp=1/.8, mu0=2.4, mu1=.6, mu2=.6, 
                          drop1=0.0001, drop2=0.0001, drop3=0.0001, fup=32/52) )  
x <- t(data.frame(res))
x <- (x<0.05)
x<- data.frame(x)
x$tally <- (1*x$p1+1*x$p2)

table(x$p1)/sims
table(x$p2)/sims
table(x$tally)/sims # neither/only 1/both
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# one data set simulation

n   <- 220
disp <- 1.3
mu0 <- 1
mu1 <- 0.65*mu0
dose <- c(rep("placebo",n),rep("trt",n))
mu   <- c(rep(mu0,n),rep(mu1,n))
drop <- c(rep(.1,n),rep(.1,n))

f <- - rexp(2*n) / log(1-drop)
length <- ifelse(f>1,1,f)

y <-  rnbinom(n*2,  p=1/(1+ mu*length* disp),      size=1/disp)  +
      rnbinom(n*2,  p=1/(1+ mu0*(1-length)*disp),  size=1/disp)

logleng  <- rep(0, n*2)

addmargins(table( y,  dose))

d <- cbind.data.frame(dose, mu, length, y, logleng)

summary(MASS::glm.nb(y~dose+offset((logleng)), data=d))



# Hilbe's neg binomial book: Table 9.5 R: Synthetic Monte Carlo negative
# no offset here!
 
# library(MASS)
# mysim <- function()
# {
#   nobs <- 220
#   x1 <- rbinom(nobs,1,prob=.5)
#   xb <- 1 + .65*x1 
#   a <- 1/1.3
#   ia <- 1/a
#   exb <- exp(xb)
#   xg <- rgamma(nobs, a, ia)
#   xbg <-exb*xg
#   nby <- rpois(nobs, xbg)
#   nbsim <-glm.nb(nby ~ x1)
#   
#   x<-summary(nbsim)
#   nbp <- x$coefficients["x1","Pr(>|z|)"]
#   
#   alpha <- nbsim$theta
#   pr <- sum(residuals(nbsim, type='pearson')^2)
#   prdisp <- pr/nbsim$df.residual
#   beta <- nbsim$coef
#   list(alpha,prdisp,nbp, beta)
# }
# B <- replicate(100, mysim())
# mean(unlist(B[1,]))
# mean(unlist(B[2,]))
# mean(unlist(B[3,])<.05)
# apply(matrix(unlist(B[4,]),2,100),1,mean)