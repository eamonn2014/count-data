

rm(list=ls())

# updating the plot in my app
 
p3 <- function(x) {formatC(x, format="f", digits=3)}
p4 <- function(x) {formatC(x, format="f", digits=4)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
p1 <- function(x) {print(formatC(x, format="f", digits=1),quote=FALSE)}
p1 <- function(x) {formatC(x, format="f", digits=1)}
p0 <- function(x) {formatC(x, format="f", digits=0)}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# try thes functions https://github.com/davidski/collector/blob/main/R/fit_distributions.R

#' Given a set of parameters describing a lognormal distribution, return
#'   the parameters of the underlying normal distribution.

lognormal_to_normal <- function(meanlog, sdlog) {
  norm_mean <- exp(meanlog + sdlog^2 / 2)
  norm_sd <- sqrt( (exp(sdlog^2) - 1) * exp(2*meanlog + sdlog^2))
  list(mean = norm_mean, sd = norm_sd)
}

#' Given parameters that describe a normal distribution, convert them back
#'   to parameters for a lognormal distribution.

#' normal_to_lognormal(normmean = 20, normsd = 3)
normal_to_lognormal <- function(normmean, normsd) {
  phi <- sqrt(normsd ^ 2 + normmean ^ 2)
  lognorm_meanlog <- log(normmean ^ 2 / phi)
  lognorm_sdlog <- sqrt(log(phi ^ 2 / normmean ^ 2))
  list(meanlog = lognorm_meanlog, sdlog = lognorm_sdlog)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mux <- 2   #  mean 
a   <- 1.5 #  standard 
 
zx  <- 0.6     # delta: hypothesised log-normal proportional change, this will be logged",
N   <- 6000    # Patients in one arm

 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  logN.2.Norm  

v <- a^2
m <- mux

phi   = sqrt(v + m^2);
mu    = log(m^2/phi)           # mean of log(Y)     
sigma = sqrt(log(phi^2/m^2))   # std dev of log(Y)  

mu;sigma

# this duplicates my code, im inpuuting lognormal 
normal_to_lognormal(normmean = mux, normsd = a)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use Normal mean and sd to get back the lognormal mean and lognormal SD!
# Norm.2.logN 

v <- sigma^2
mu <- mux

logN.mu <- exp(mu+0.5* v)
logN.SD <- logN.mu*sqrt(exp(v)-1)

logN.mu;logN.SD

# this duplicates my code
lognormal_to_normal(meanlog=mu, sdlog=sigma)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate data

sd    <- sigma
mu    <- mu
delta <- log(zx)

A <- lapply(1:N, function(i) rlnorm(1, meanlog=mu,       sdlog=sd))
B <- lapply(1:N, function(i) rlnorm(1, meanlog=mu+delta, sdlog=sd))

A <- unlist(A)
B <- unlist(B)

A.GM <- exp(mean(log(A)))
B.GM <- exp(mean(log(B)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow=c(2,2))

x<-seq(-8+mu,+8+mu,by=0.02)

br <- ifelse(N<100, N, 100) # breaks in histogram

#~~~~~~~~~~~~~~~~~~~
MASS::truehist(A, 
              yaxt='n' ,
              nbins=br,  axes=FALSE,
               main=paste0("N=",N," realisations from log-normal mu=",mux,", SD=",p3(sigma),
                           ",\n Geo.Mean=",p3(A.GM)), col = "#75AADB", 
              border = "white", xlab="Original scale")

curve(dlnorm(x, meanlog = mu, sdlog =sigma, log = FALSE), add=TRUE)
Axis(side=1, labels=TRUE)
Axis(side=2, labels=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MASS::truehist(B, 
                yaxt='n' ,
                nbins=br,  axes=FALSE,
               main=paste0("N=",N," realisations true log-normal mu=",mux*zx,", SD=", p3(sigma),
                           ", \nGeo.Mean=",p3(B.GM) ,"\nTreatment effect ",mux*zx," / ",mux,"=",zx), 
               col = "red", border = "white", xlab="Original scale")
 
curve(dlnorm(x, meanlog = mu+delta, sdlog =sigma, log = FALSE), add=TRUE)
 Axis(side=1, labels=TRUE)
Axis(side=2, labels=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MASS::truehist(log(A), yaxt='n' , nbins=br,  axes=FALSE,
               main=paste0("N=",N," realisations true Normal \nmu=",p3(mu),", SD=",
                           p3(sd),""), col = "#75AADB", border = "white",  xlab="log scale")
curve(dnorm(x,mu,sd), add=TRUE)
Axis(side=1, labels=TRUE)
Axis(side=2, labels=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MASS::truehist(log(B), yaxt='n' , nbins=br, axes=FALSE, 
               main=paste0("N=",N," realisations true Normal mu=",p3(mu+delta),
                           ", SD=",p3(sd)," \nTreatment effect ", p2(mu+delta)," - ",p3(mu),"=",
                           p3(log(zx))), col = "red", border = "white", xlab="log scale" )
curve(dnorm(x,mu+delta,sd), add=TRUE)
Axis(side=1, labels=TRUE)
Axis(side=2, labels=FALSE)

par(mfrow=c(1,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# boxplots

conf    <- 0.9
qclevel <- 0.85

d <- data.frame(x =c(rep( "A" , N), rep( "B" , N)), 
                y=c(A,B))

res <- t.test(log(d$y)~d$x, var.equal = TRUE, conf.level=conf) # t test

ci <- exp(-res$conf.int)
hi <- round(ci[1],3)
lo <- round(ci[2],3)

qc <- ifelse(((qclevel<=lo) && (hi<=1/qclevel)), paste0("PASS spec. of ", qclevel), paste0("FAIL spec. of ", qclevel))
GmeanR <- round(exp(res$estimate[1]),3) # geometric means
GmeanT <- round(exp(res$estimate[2]),3) # geometric means

# %CV get the sd of the logged data and exponenitate
cv1 <-exp(tapply(log(d$y), d$x, sd))
GmeanRcv <-  round((((cv1[1])-1)*100),3)
GmeanTcv <-  round((((cv1[2])-1)*100),3)

ratio. <- exp( (res$estimate)[2] - (res$estimate)[1])
ratio. <- round(ratio.,3)
require(ggplot2)
 

# run a ggplot only for the purposes of get the ticks for the y axes
d$Y <-  (d$y)  # choose response

gx <- ggplot(d, aes(x=x, y=Y)) +
  geom_point(aes(fill=x), size=3, shape=21, colour="grey20",
             position =position_jitter(width=.22, height=0))+
  geom_boxplot(outlier.colour=NA, fill=NA,colour="grey20") 
  
  
print(gx)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### log the response and use exponentiated ticks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # now use the tick info on ggplot of log response
   d$Y <-  log(d$y)  # choose response
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # try and control the y axis
   if(min(exp(d$Y) ) < 0.01) {
   ticks= log(c(0.001, 0.01, .1 ,1, 10, 100))
   } else {
   ticks= log(c(0.01, .1 ,1, 10, 100))
   }
   
   if( min(exp(d$Y)) > 0.1) {
     ticks= log(c(0.1, .1 ,1, 10, 100))
   } else {
     ticks= log(c(0.01, .1 ,1, 10, 100))
   }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   labs <- exp(ticks)  # exponentiate the labels
   
   g1 <- ggplot(d, aes(x=x, y=Y)) +
     geom_point(aes(fill=x), size=3, shape=21, colour="grey20",
                position =position_jitter(width=.22, height=0))+
     geom_boxplot(outlier.colour=NA, fill=NA,colour="grey20") +
      
    scale_y_continuous(
      
     limits=c(   ticks[1],  ticks[length(ticks)]  ), 
        breaks= ticks,    # this is where the values go
        labels= labs) +   # these are labels
     
     labs(title=paste0("QC=", qc)) +
     labs(subtitle = 
            paste0("Geometric mean A=" ,GmeanR,"; %CV=",GmeanRcv,
                   "\nGeometric mean B=" ,GmeanT,"; %CV=",GmeanTcv,
                   "\nRatio B/A " ,ratio.,", ",conf*100,"%CI (",lo,", ",hi,")")) +
     
     annotate("text", x=0.5, y=max(ticks), label=paste0("n=", table(d$x)[1][[1]])) +
     annotate("text", x=1.5, y=max(ticks), label=paste0("n=", table(d$x)[2][[1]]))          
   
   print(g1)
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   







