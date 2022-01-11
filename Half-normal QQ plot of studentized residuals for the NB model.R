
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Atkinson (1981, 1987) suggested a more robust and useful version of these QQ plots: half
# normal plots, with simulated confidence envelopes. Friendly page 481


    rm(list=ls())
    require(tidyverse)
    require(MASS)
    library(car)
    library(vcd)
     
    
    # pop parameters
    n  <- 1000
    k  <- 1.3   # this is k, alpha=1/k
    1/k         # alpha reported as theta in neg binomial
    mu0 <- 1
    mu1 <- mu0*0.75
    
    fup   <- 1
    drop1 <- 0.1
    drop2 <- 0.1
     
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dose <- c(rep("placebo",n),rep("trt",n))
    mu   <- c(rep(mu0,n),rep(mu1,n))
    drop <- c(rep(.1,n),rep(.1,n))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # discontinuations and follow up
    f        <- - rexp(2*n) / log(1-drop)
    length   <- ifelse(f>1,1,f)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # rnbinom approach
    y <-  rnbinom(n*2,  p=1/(1+ mu*    length* k),  size=1/k)  +
          rnbinom(n*2,  p=1/(1+ mu0*(1-length)*k),  size=1/k)
    
    logleng  <- rep(0, n*2)
    
    summary(f<-glm.nb(y~dose+offset(logleng) ))



  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  par(mfrow=c(1,2))
  qqPlot(rstudent(f), xlab="Normal quantiles", ylab="Studentized residuals")
  
  
  # michael friendly
  # Figure 11.38: Half-normal QQ plot of studentized residuals for the NB model fit to the PhdPubs
  # data. The reference line and confidence envelope reflect the mean and (2.5%, 97.5%) quantiles of
  # {fig:phd-halfnorm} the simulation distribution under the negative-binomial model for the same data
  
  # First, calculate the sorted absolute values of the residuals |r|(i) and their 
  # expected normal values, |z|(i) . The basic plot will be plot(expected, observed).
  
  
  observed <- sort(abs(rstudent(f)))
  n <- length(observed)
  expected <- qnorm((1:n + n - 1/8)/(2*n + 1/2))
  
  #
  # Then, use simulate() to generate S = 100 simulated response vectors 
  # around the fitted values in the model. Here this uses the negative-binomial
  # random number generator (rnegbin()) with the same dispersion value 
  # estimated in the model. The result, called sims here, 
  # is a data frame of n rows and S = 100 columns, named sim_1, sim_2, ....
  
  S <- 100
  sims <- simulate(f, nsim=S)
 
  # The next step is computationally intensive, because we have to fit the 
  # NB model S = 100 times and a little bit tricky, because we need to use
  # the same model formula as the original, but with the simulated y. We 
  # first define a function resids to do this for a given y, and then use
  # apply to calculate them all. To save computing time,
  # the coefficients from the original model are used as starting values.
  
  resids <- function(y)
    rstudent(glm.nb(y ~ dose + offset(logleng),
                    start=coef(f)))
  
  # fit models to each simulated y and collect student residuals
  z <- apply(sims,2,resids)
  z <- abs(z)            # absolute resids
  z <- apply(z, 2, sort) # sort each column
  
  # now calculate the mean and quantiles 
  envelope <- 0.95
  mean <- apply(z, 1, mean)
  lower <- apply(z, 1, quantile, prob=(1 - envelope)/2)
  upper <- apply(z, 1, quantile, prob=(1 + envelope)/2)
  
  
  plot(expected, observed,
       xlab='Expected value of half-normal order statistic',
       ylab='Absolute value of studentized residual')
  lines(expected, mean, lty=1, lwd=2, col="blue")
  lines(expected, lower, lty=2, lwd=2, col="red")
  lines(expected, upper, lty=2, lwd=2, col="red")
  # identify(expected, observed, labels=names(observed), n=3)
   
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  res <- rstudent(f)
  plot(density(res), lwd=2, col="blue",
         main="Density of studentized residuals")
   rug(res)
  
    #why the bimodality?
   plot(jitter(log(y+1), factor=1.5), res,
            xlab="log (y+1)", ylab="Studentized residual")
  
   ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  