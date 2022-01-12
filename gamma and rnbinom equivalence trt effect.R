
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
n  <- 100000
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
 
f        <- - rexp(2*n) / log(1-drop/fup) # scale according to follow up!
length   <- ifelse(f > fup, fup, f)  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# rnbinom approach
y <-  rnbinom(n*2,  p=1/(1+ mu*    length* k),  size=1/k)  +
      rnbinom(n*2,  p=1/(1+ mu0*(1-length)*k),  size=1/k)

logleng  <- rep(log(fup), n*2) 
 
summary(f<-glm.nb(y~dose+offset(logleng) ))


#~assumption checking
# https://stats.stackexchange.com/questions/70558/diagnostic-plots-for-count-regression
performance::check_model(f)
library(vcd)


Ord_plot(y)
distplot(y, type="poisson")
distplot(y, type="nbinom")

library(AER)
deviance(f)/f$df.residual
dispersiontest(f)

library(car)
influencePlot(f)

res <- residuals(f, type="deviance")
plot(log(predict(f)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Simulate the gamma-shaped subject effect, same for each arm
s <-  rgamma(n, shape = 1/k, scale = k)
s1 <- rgamma(n, shape = 1/k, scale = k)

lambda1 <- mu   * s  ## Simulate rate
lambda2 <- mu0 * s1  ## Simulate rate

y_0 <- rpois(n = n*2, lambda = (length*lambda1))  +
       rpois(n = n*2, lambda = ((1-length)*lambda2))  
 
logleng  <- rep(log(fup), n*2)  # all patients same f up

summary(f1<-glm.nb(y_0~dose+offset(logleng) ))


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






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~simulation
#https://www.jamesuanhoro.com/post/2018/05/07/simulating-data-from-regression-models/

# column is a single simulation
nn <-100
sim.default <- simulate(f, nn)
dim(sim.default)

z <- sample(1:nn,1)
ys <- sim.default[,z]
summary(glm.nb(ys~dose+offset(logleng) ))


 

#https://data.library.virginia.edu/simulating-data-for-count-models/
require(countreg)
countreg::rootogram(f)

countreg::rootogram(f, style = "standing",   main = "Standing")
countreg::rootogram(f, style = "hanging",     main = "Hanging")
countreg::rootogram(f, style = "suspended",   main = "Suspended")
 
rootogram(f, style = "hanging",     main = "Hanging", scale="raw")
 
## inspect output (without plotting)
r <- rootogram(f, plot = FALSE)
r









coefs <- mvrnorm(n = 10000, mu = coefficients(f), Sigma = vcov(f))

# coef
coefficients(f)
colMeans(coefs) #
#SEs
sqrt(diag(vcov(f)))
apply(coefs, 2, sd) # 



# # One row per case, one column per simulated set of coefficients
# sim.dat <- matrix(nrow = n, ncol = nrow(coefs))
# fit.p.mat <- model.matrix(f) # Obtain model matrix
# # Cross product of model matrix by coefficients, exponentiate result,
# # then use to simulate Poisson-distributed outcome
# for (i in 1:nrow(coefs)) {
#   sim.dat[, i] <- rpois(n, exp(fit.p.mat %*% coefs[i, ]))
# }
# rm(i, fit.p.mat) # Clean house


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~trying predict and simulate

fp <- predict(f , newdata=data.frame(dose=dose))


# not sure if this is right


  # set.seed(4)
  n <- 100
  
  dose <- c(rep("placebo",n),rep("trt",n))
  
  k <- f$theta
  intercept <- as.vector(f$coefficients)[1]
  beta      <- as.vector(f$coefficients)[2]
  
  #p.out <- predict(f, type = "response")
  
  s <-  rgamma(n, shape = 1/k, scale = k)
  s1 <- rgamma(n, shape = 1/k, scale = k)
  
  lambda1 <- exp(intercept)   * s  ## Simulate rate
  lambda2 <- exp(intercept + beta) * s1  ## Simulate rate
  
  fx        <- - rexp(2*n) / log(1-drop)
  length   <- ifelse(fx>1 ,1, fx)
  
  
  y_0 <- rpois(n = n, lambda = (length*    lambda1))  +
         rpois(n = n, lambda = ((1-length)*lambda2))  
  
  y_1 <- rpois(n = n, lambda = (length*    lambda2))  +
         rpois(n = n, lambda = ((1-length)*lambda2))  
  
  dat2 <- data.frame(counts =c(y_0, y_1), dose=dose)
   
  ggplot(dat2, aes(x = counts)) +
    geom_bar() +
   facet_wrap(~dose)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  

  
###
# library(MASS)
# mysim <- function()
# {
#   nobs <- 500
#   x1 <-runif(nobs)
#   x2 <-runif(nobs)
#   xb <- 2 + .75*x1 - 1.25*x2
#   a <- .5
#   ia <- 1/.5
#   exb <- exp(xb)
#   #xg <- rgamma(nobs, a, ia) # intercept is wrong with this code
#   xg  <- rgamma(nobs, shape = a, scale = ia)
#   xbg <-exb*xg
#   nby <- rpois(nobs, xbg)
#   nbsim <-glm.nb(nby ~ x1 + x2)
#   alpha <- nbsim$theta
#   pr <- sum(residuals(nbsim, type="pearson")^2)
#   prdisp <- pr/nbsim$df.residual
#   beta <- nbsim$coef
#   list(alpha,prdisp,beta)
# }
# B <- replicate(10, mysim())
# mean(unlist(B[1,])) #theta , this is alpha
# mean(unlist(B[2,])) #Pearson-dispersion statistic
# apply(matrix(unlist(B[3,]),3,100),1,mean) # intercept x1 x2
# 





