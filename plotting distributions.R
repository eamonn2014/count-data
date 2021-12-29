

# 
# 
# x <- sort(rnorm(100,0,0.5))
# h <- hist(x,plot=FALSE)
# dens1 <-  h$counts/sum(h$counts)
# dens2 <- dnorm(x,0,0.5)
# 
# hist(x,probability=TRUE,breaks="fd",ylim=c(0,1))
# lines(h$mids,dens1,col="red")
# lines(x,dens2,col="darkgreen")
 
 
  
  N=500
  p=runif(1,0,1)
  #tmp <- sample(c(runif(1,0,1),c(runif(1,1,1000))),1)
  s <- runif(1,0,4)   # hilbe stated alpha not greater than 4 usually
  #s=round(tmp,2)
  mu = s/(p) - s
  
   x <-0:round(mu*2,0)
   A <- rnbinom(N, size = s, prob = p)                          # Draw N neg binomially distributed values
   m <- max(c(dpois(x, lambda=mu),dnbinom(x, size=s, mu=mu)))   # help with plotting 
   
  hist(A, prob=TRUE, breaks=100, ylim=c(0,m),
       main=paste0("Negative binomial N=",N,", prob=",round(p,3),", alpha=",round(s,3),", mu=",round(mu,3) ) )
  curve(dnorm(x, mu, sqrt(mu)), lwd=2, col="red",  add=T) # 'x' mandatory arg
  
  points(x, dpois(x, lambda=mu), pch=19, col="darkgreen")   # poisson
  
  points(x, dnbinom(x, size=s, mu=mu), pch=19, col="purple")  # fitted neg binomial
  
  legend(x = "topright",                           # Position
         legend = c("norm", "poisson","neg bin"),  # Legend texts
         lty = c(1),                               # Line types
         col = c("red","darkgreen","purple"),      # Line colors
         lwd = 2,title.adj = 0.95,
             box.lty = 0)                          # Line width
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  N=500
  p=runif(1,0,1)
  s <- runif(1,0,4)   # hilbe stated alpha not greater than 4 usually
  mu = s/(p) - s
   
  A <- rnbinom(N, size = s, prob = p)          # Draw N neg binomially distributed values
  x <- as.numeric(names(table(A)))             # use for plotting
  m <- max(c(dpois(x, lambda=mu),dnbinom(x, size=s, mu=mu))) # use for plotting
  
  plot(dnorm(x, mu, sqrt(mu)), lwd=2, col="red",type='l', lty=1, xlim=c( 0, mu+8*sqrt(mu)), ylim=c(0,m),
  
  main=paste0("Negative binomial N=",N,", prob=",round(p,3),", alpha=",round(s,3),", mu=",round(mu,3) ) )
  
  lines(x, dpois(x, lambda=mu), pch=19, col="darkgreen")
  
  lines(x, dnbinom(x, size=s, mu=mu), pch=19, col="purple")
  
  legend(x = "topright",                           # Position
         legend = c("norm", "poisson","neg bin"),  # Legend texts
         lty = c(1),                            # Line types
         col = c("red","darkgreen","purple"),      # Line colors
         lwd = 2,title.adj = 0.95,
         box.lty = 0)                          # Line width
  
  
  
  
  
  
  
  
  
  
  
  
 #  
 #  lines(A, col="red")
 #  
 # #https://math.stackexchange.com/questions/2412983/plotting-in-r-probability-mass-function-for-a-poisson-distribution
 #  
 #  x = 0:20;  pdf = dpois(x, 6)
 #  plot(x, pdf, type="h", lwd=3, col="blue", 
 #       main="PDF of POIS(6) with Approximating Normal Density")
 #  abline(h=0, col="green2")
 #  curve(dnorm(x, 6, sqrt(6)), lwd=2, col="red", add=T) # 'x' mandatory arg