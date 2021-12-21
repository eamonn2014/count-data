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
  
  # analyse with neg. binomial model
  x <- summary(MASS::glm.nb(y~dose+offset((logleng))))
  
  # collect p-values
  p <-  x$coefficients["dosetrt","Pr(>|z|)"]
  
  return(p)
  
}

res <- replicate(1000, nb.power())  #.77
mean(res<0.05)
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end make a function
res <- replicate(1000, nb.power(mu0=1.5,n=120))  #.59
mean(res<0.05)
res <- replicate(1000, nb.power(mu0=.9,n=120))   #.49
mean(res<0.05)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# one data set simulation

# n   <- 2200
# disp <- 1.3
# mu0 <- 1
# mu1 <- 0.65*mu0
# dose <- c(rep("placebo",n),rep("trt",n))
# mu   <- c(rep(mu0,n),rep(mu1,n))
# drop <- c(rep(.1,n),rep(.1,n))
# 
# f <- - rexp(2*n) / log(1-drop)
# length <- ifelse(f>1,1,f)
# 
# y <-  rnbinom(n*2,  p=1/(1+ mu*length* disp),      size=1/disp)  +
#   rnbinom(n*2,  p=1/(1+ mu0*(1-length)*disp),  size=1/disp)
# 
# logleng  <- rep(0, n*2)
# 
# addmargins(table( y,  dose))
# 
# d <- cbind.data.frame(dose, mu, length, y, logleng)
# 
# summary(MASS::glm.nb(y~dose+offset((logleng)), data=d))
