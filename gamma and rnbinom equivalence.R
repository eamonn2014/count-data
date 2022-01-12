
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# gamma and rnbinom equivalence 

  rm(list=ls())
  require(tidyverse)
  require(MASS)
  
  roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # pop parameters
  n  <- 1000
  k  <- 1.3   # this is k, alpha=1/k
  1/k         # alpha reported as theta in neg binomial
  mu <- 2
  
  fup  <- 2
  drop <- 0.4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # discontinuations and follow up
   
  f        <- - rexp(n) / log(1-drop/fup) # scale according to follow up!
  length   <- ifelse(f > fup, fup, f)   
  logleng  <- log(length)  # offset
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  # rnbinom approach
  mu   <-   c(rep(mu,n) )
  y <-  rnbinom(n,  p=1/(1+ mu*length* k),  size=1/k)  
  summary(glm.nb(y~1+offset(logleng)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Simulate the gamma-shaped subject effect
  s <- rgamma(n, shape = 1/k, scale = k)
  lambda1 <- mu * s  ## Simulate rate
  y_0 <- rpois(n = n, lambda = (length*lambda1))  # Simulate counts, added length here
  summary(glm.nb(y_0~1 +offset(logleng)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # plot
  
  trt <- y
  pla <- y_0
  
  # plot counts
  par(mfrow=c(2,2))
  u <- roundUpNice(max(trt %>% table, pla %>% table))
  trt %>% table %>% barplot(ylim=c(0,u)) 
  pla %>% table %>% barplot(ylim=c(0,u)) 
  trt %>% table %>% prop.table
  pla %>% table %>% prop.table
  # par(mfrow=c(1,1))
  
  # plot proportion
  # par(mfrow=c(1,2))
  data_perc <- t(prop.table(table(trt)))  * 100     # Convert data to probability table
  data_perc1 <- t(prop.table(table(pla))) * 100     # Convert data to probability table
  u <- roundUpNice(max(data_perc,data_perc1))
  barplot(data_perc,  ylab = "Percent", main ="rnbinom", ylim=c(0,u))
  barplot(data_perc1, ylab = "Percent", main ="gamma",   ylim=c(0,u))
  par(mfrow=c(1,1))
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  
   
  
  
  
  
  
  
  