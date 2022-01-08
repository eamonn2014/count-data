
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# gamma and rnbinom equivalence 

rm(list=ls())
require(tidyverse)
require(MASS)

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# pop parameters
n  <- 10000
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
 
summary(MASS::glm.nb(y~dose+offset(logleng) ))
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Simulate the gamma-shaped subject effect, same for each arm
s <-  rgamma(n, shape = 1/k, scale = k)
s1 <- rgamma(n, shape = 1/k, scale = k)

lambda1 <- mu   * s  ## Simulate rate
lambda2 <- mu0 * s1  ## Simulate rate

y_0 <- rpois(n = n*2, lambda = (length*lambda1))  +
       rpois(n = n*2, lambda = ((1-length)*lambda2))  
 
       logleng  <- rep(0, n*2)  # all patients have same length follow up

summary(glm.nb(y_0~dose+offset(logleng) ))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trt <- y[dose %in% "trt"]
pla <- y[!dose %in% "trt"]
# examine
par(mfrow=c(2,2))
u <- roundUpNice(max(trt %>% table, pla %>% table))
trt %>% table %>% barplot(ylim=c(0,u)) 
pla %>% table %>% barplot(ylim=c(0,u)) 
trt %>% table %>% prop.table
pla %>% table %>% prop.table
#par(mfrow=c(1,1))

#par(mfrow=c(1,2))
data_perc <- t(prop.table(table(trt)))  * 100     # Convert data to probability table
data_perc1 <- t(prop.table(table(pla))) * 100     # Convert data to probability table
u <- roundUpNice(max(data_perc,data_perc1))
barplot(data_perc,  ylab = "Percent", main ="rnbinom", ylim=c(0,u))
barplot(data_perc1, ylab = "Percent", main ="gamma",   ylim=c(0,u))   
par(mfrow=c(1,1)) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trt <- y_0[dose %in% "trt"]
pla <- y_0[!dose %in% "trt"]
# examine
par(mfrow=c(2,2))
u <- roundUpNice(max(trt %>% table, pla %>% table))
trt %>% table %>% barplot(ylim=c(0,u)) 
pla %>% table %>% barplot(ylim=c(0,u)) 
trt %>% table %>% prop.table
pla %>% table %>% prop.table
#par(mfrow=c(1,1))

#par(mfrow=c(1,2))
data_perc <- t(prop.table(table(trt)))  * 100     # Convert data to probability table
data_perc1 <- t(prop.table(table(pla))) * 100     # Convert data to probability table
u <- roundUpNice(max(data_perc,data_perc1))
barplot(data_perc,  ylab = "Percent", main ="rnbinom", ylim=c(0,u))
barplot(data_perc1, ylab = "Percent", main ="gamma",   ylim=c(0,u))   
par(mfrow=c(1,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~












