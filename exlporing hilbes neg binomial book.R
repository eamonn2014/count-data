
#  negative binomial heterogeneity (ancillary) parameter 
#  The only restriction placed on α is that it take positive rational values, rarely above 4 # hilbe p 190
#  there are some error in hilbe code that  have corrected - rgamma specificiation
require(MASS) # p225
nb2_syn <- function(nobs = 100, off = 0,
                    alpha = .5,
                    xv = c(2, 0.75, -1.5)) {
  p <- length(xv) - 1
  X <- cbind(1, matrix(rnorm(nobs * p), ncol = p))
  xb <- X %*% xv
  a <- alpha
  ia <- 1/a
  exb <- exp(xb + off)
#  xg <- rgamma(nobs,a,  a, ia)   # this is wrong!!
  xg <- rgamma(nobs, shape = a, scale = ia) 
  xbg <-exb*xg
  nby <- rpois(nobs, xbg)
  out <- data.frame(cbind(nby, X[,-1]))
  names(out) <- c("nby", paste("x", 1:p, sep=""))
  return(out)
   
}
# --------------------------------------------------
sim.data <- nb2_syn(nobs = 50000, alpha=.5, xv = c(2, .75,
                                                   -1.25))
mynb <- glm.nb(nby ~ ., data = sim.data)

f <- summary(mynb)
f
# i <- exp(f$coeff)[1.1]
#  1/exp(-f$coeff)[1.1]
# nd = data.frame(x1=0,x2=0)
# nd = data.frame(x1=mean(sim.data$x1),x2=mean(sim.data$x2))
# predict(mynb,newdata=nd,type="response")


 
 require(MASS)
 # choose one or other sample size
 n <- 20000
 # n <- 10000
 # generating values
 beta <- c(0.5,0.8,-0.4) # regression parameters
 sigma <- 0.7 # overdispersion parameter
 # sample 100 datasets for given n
 beta.sample.data <- matrix(NA,100,3)
 X1 <- rnorm(n)
 X2 <- rnorm(n)
 #X2 <- c(rep(0,n/2),rep(1,n/2))
 mu <- exp(beta[1]+beta[2]*X1+beta[3]*X2) # regression means
 y <- rnbinom(n=n,mu=mu,size=1/sigma) # generate responses
 NB <- glm.nb(y ~X1+X2)
 NB
 


# Hilbe p229, takes a few  minutes

library(MASS)
mysim <- function()
{
  nobs <- 50000
  x1 <-runif(nobs)
  x2 <-runif(nobs)
  xb <- 2 + .75*x1 - 1.25*x2
  a <- .5
  ia <- 1/.5
  exb <- exp(xb)
  #xg <- rgamma(nobs, a, ia) # intercept is wrong with this code
  xg  <- rgamma(nobs, shape = a, scale = ia)
  xbg <-exb*xg
  nby <- rpois(nobs, xbg)
  nbsim <-glm.nb(nby ~ x1 + x2)
  alpha <- nbsim$theta
  pr <- sum(residuals(nbsim, type="pearson")^2)
  prdisp <- pr/nbsim$df.residual
  beta <- nbsim$coef
  list(alpha,prdisp,beta)
}
B <- replicate(10, mysim())
mean(unlist(B[1,])) #theta , this is alpha
mean(unlist(B[2,])) #Pearson-dispersion statistic
apply(matrix(unlist(B[3,]),3,100),1,mean) # intercept x1 x2