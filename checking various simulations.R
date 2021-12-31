# this is useful
# correlated poisson
# correlated neg binomial 
# two methods of correlation

# hilbe page 186
# Note that the geometric model is a variety of negative binomial, i.e. a negative
# binomial with α = 1. The Poisson may also be included as a negative binomial,
# i.e. α = 0. 

# intercepts is off but the beta is good with high dispersion (<1)!

rm(list=ls())

library(MASS)
library(lme4)

# parameters
n <- 220
r <- .75  
mu0 <- 10 
mu1 <- .75 # parameters 
disp <- .8 # use this for neg binomial - big dispersion means we will approx poisson

L1 <- mu0
L2 <- mu0*mu1     # change by the intervention

# code from business paper

GenerateMultivariatePoisson<-function(p, samples, R, lambda){
  normal_mu=rep(0, p)
  normal = mvrnorm(samples, normal_mu, R)
  unif=pnorm(normal)
  pois=t(qpois(t(unif), lambda))
  return(pois)
}

# execute and generated correlated poisson
d <-  GenerateMultivariatePoisson(p=2, samples=n, R=matrix(c(1, r, r, 1), 2, 2), lambda=c(L1,L2))


mean(d[,1]) # check mean rates
mean(d[,2])
cor(d)

 # create a long dataset
my_data <- data.frame( 
  group = rep(c("before", "after"), each = n),
  counts = c(d[,1],  d[,2]),
  ID=rep(1:n,2)
)


# examine
require(tidyverse)
par(mfrow=c(1,2))
d[,1] %>% table %>% barplot() #quick and dirty
d[,1] %>% table %>% prop.table
d[,2] %>% table %>% barplot() #quick and dirty
d[,2] %>% table %>% prop.table
par(mfrow=c(1,1))


# better plot
par(mfrow=c(1,2))

z<- d[,1]
breakz= 15  # why 15?
h<-hist(z,breaks=breakz, plot=F)
xhist<-c(min(h$breaks),h$breaks)
yhist<-c(0,h$density,0)
xfit<-seq(min(z),max(z))

yfit<-dpois(xfit,lambda=L1)  # mu=mean(z)?
plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)),
     main=paste0("Poisson pdf and hist N=",n,", \nrate=",round(L1,3),"" ) )
lines(xfit,yfit, col='red')


z<- d[,2]
breakz= 15  # why 15?
h<-hist(z,breaks=breakz, plot=F)
xhist<-c(min(h$breaks),h$breaks)
yhist<-c(0,h$density,0)
xfit<-seq(min(z),max(z))

yfit<-dpois(xfit,lambda=L2)  # mu=mean(z)?
plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)),
     main=paste0("Poisson pdf and hist N=",n,", \nrate=",round(L2,3),"" ) )
lines(xfit,yfit, col='red')

par(mfrow=c(1,1))


 
# analyse
A <- glmer(counts ~ group + (1|ID), data=my_data, family="poisson")
x <- summary(A)
(rate.reduction <- 1-1/exp(x$coefficients[2,1]))
exp(x$coefficients)


####################################################################################################

# based on mu size (dispersion) create functions to swap between prob and mu

mu2prb <- function(mu, size){
  
  prob = size/(size+mu) 
  print(prob)
}

prb2mu <- function(prb, size){
  
  mu = size/(prb) - size 
  
  print(mu) 
}



p1 <- mu2prb(10,disp)   # dont think there is any need for this as we could prob use mu
p2 <- mu2prb(7.5,disp)

library(SimCorrMix)  # try this library for generating correlated ned binomial
# seed will default so I change each time
s <- corrvar2(n =n, 
              k_nb = 2,  
              size = c(disp, disp), prob = c(p1,p2), 
              rho = matrix(c(1, 0.75, 0.75, 1), 2, 2), p_zinb = 0, seed=sample(1000:9999,1)) # 

apply(s$Y_n,2,mean)  # check mean rate

m1 <- prb2mu(p1, disp)  # rate 1
m2 <- prb2mu(p2, disp)  # rate 2

cor(s$Y_n[,1:2]) # check corr

# examine
par(mfrow=c(1,2))
s$Y_n[,1] %>% table %>% barplot() #quick and dirty
s$Y_n[,1] %>% table %>% prop.table
s$Y_n[,2] %>% table %>% barplot() #quick and dirty
s$Y_n[,2] %>% table %>% prop.table
par(mfrow=c(1,1))

# we are here trying ot plot better
# 
# hist(s$Y_n[,1] , prob=TRUE, breaks=220, #ylim=c(0,1),
#      main=paste0("Negative binomial N=",n,", prob=",round(p1,3),", disp=",round(disp,3),", mu=",round(m1,3) ) )
# #curve(dnorm(x, mu, sqrt(mu)), lwd=2, col="red",  add=T) # 'x' mandatory arg
# #
# # barplot(table(s$Y_n[,1]))
# points(dnbinom(x=2:20, size=disp, mu=m1), lwd=2, col="red",  add=T) # 'x' man

# better examination

# https://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf p15
# A graphical technique to evaluate the goodness of fit can be drawing pdf curve and histogram together
par(mfrow=c(1,2))

  z<- s$Y_n[,1]
  breakz= 15  # why 15?
  h<-hist(z,breaks=breakz, plot=F)
  xhist<-c(min(h$breaks),h$breaks)
  yhist<-c(0,h$density,0)
  xfit<-seq(min(z),max(z))
  
  yfit<-dnbinom(xfit,mu=m1, size=disp)  # mu=mean(z)?
  plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)),
       main=paste0("Neg binomial pdf and hist N=",n,", \nprob=",round(p1,3),", dispersion=",round(disp,3),", mu=",round(m1,3) ) )
  lines(xfit,yfit, col='red')
  
  
  z<- s$Y_n[,2]
  h<-hist(z,breaks=breakz, plot=F)
  xhist<-c(min(h$breaks),h$breaks)
  yhist<-c(0,h$density,0)
  xfit<-seq(min(z),max(z))
  yfit<-dnbinom(xfit,mu=m2, size=disp)
  plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)), 
  main=paste0("N=",n,", prob=",round(p2,3),", \ndispersion=",round(disp,3),", mu=",round(m2,3) ) )
  lines(xfit,yfit, col='red')

par(mfrow=c(1,1))


# create long data frame
my_data <- data.frame( 
  group = rep(c("first", "second"), each = n),
  counts = c(s$Y_n[,1] ,  s$Y_n[,2] ),
  ID=rep(1:n, 2)
)

# negative binomial and possion regression, do we recover parameters?
library(lme4)
A <- NULL
A <- glmer.nb(counts ~ group + (1|ID), data=my_data ) # slow?
x <- summary(A)
# getME(A, "glmer.nb.theta")                            # what is this ? thought it would be dispersion
exp(x$coefficients)

x <- summary( f<- glmer(counts ~ group + (1|ID), data=my_data, family="poisson"))
exp(x$coefficients)

anova(A,f)

# extra
# ben bolker https://stackoverflow.com/questions/22842017/model-checking-and-test-of-overdispersion-for-glmer
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(f)  # check poisson model
overdisp_fun(A)  # check nb model


# hilbe neg binomial book p175

poi <- f
mu <-predict(poi, type="response")
z <- ((my_data$counts - mu)^2 - my_data$counts)/ (mu * sqrt(2))
zscore <- lm(z ~ 1)
summary(zscore)

# hilbe p181
library(MASS)
xb <- 3
exb <- exp(xb)
yp <-rpois(50000, exp(xb))
p3 <-glm(yp ~ 1, family=poisson)
summary(p3)
exb
 
a <- .33333333  # alpha
ia <- 1/.33333333  # dispersion
xg <- rgamma(50000, shape=a,   scale=ia)
xbg <-exb*xg
ynb <- rpois(50000, xbg)
nb3 <-glm.nb(ynb ~ 1)
summary(nb3)













