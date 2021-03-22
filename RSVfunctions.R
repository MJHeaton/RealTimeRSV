###
## Functions for the RSV analysis
###

xi.to.delta <- function(xis){
  return(plogis(xis))
}

delta.to.t <- function(deltas){
  dseq <- 1:length(deltas)
  cumsum(sapply(dseq,function(x){24*deltas[x]*prod(1-deltas[dseq<x])}))
}

rsv_curve <- function(t1,t2,t3,t4,time, beta0, beta1){
  
  out <- rep(beta0,length(time))
  out[t1<time & time <=t2] <- beta0 + beta1*(time[t1<time & time <=t2]-t1)
  out[t2<time & time <=t3] <- beta0 +beta1*(t2-t1)
  m <- -(beta1*(t2-t1))/(t4-t3)
  b <- beta0-m*t4
  out[t3<time & time <=t4] <- m*time[t3<time & time <=t4]+b
  return(out)
  
}

mult_rsv_curves <- function(t1,t2,t3,t4,time,beta0,beta1){
  nyears <- length(beta1)
  return(c(sapply(1:nyears,function(x){rsv_curve(t1[x],t2[x],t3[x],t4[x],time,beta0,beta1[x])})))
}

curve.to.prob <- function(curve.val){
  return(plogis(curve.val))
}

####### Test the functions ############
library(magrittr) #for piping %>%
OC <- read.csv("theOC.csv")[,-1]
OC <- OC[which(OC$Strt=="2003-07-01"):which(OC$End=="2013-06-30"),]
OC <- within(OC,{
  cycle <- rep(1:(nrow(OC)/24),each=24)
})
num.years <- max(OC$cycle)
xi <- matrix(rnorm(num.years*4,-1,.5),nrow=num.years,ncol=4) #Random time periods
delta <- apply(xi,1,xi.to.delta) %>% t() #convert to delta
ts <- apply(delta,1,delta.to.t) %>% t() #convert to t's
t1 <- ts[,1]
t2 <- ts[,2]
t3 <- ts[,3]
t4 <- ts[,4]
time.seq <- 1:24
the.year <- 1
beta0 <- -3
beta1 <- exp(rnorm(num.years,log(.1),.1))
plot(time.seq,
     rsv_curve(t1[the.year],t2[the.year],t3[the.year],t4[the.year],time.seq,beta0,beta1[the.year]),
     type="l")
plot(time.seq,
     curve.to.prob(rsv_curve(t1[the.year],t2[the.year],t3[the.year],t4[the.year],time.seq,beta0,beta1[the.year])),
     type="l")
plot(mult_rsv_curves(t1,t2,t3,t4,time.seq,beta0,beta1),type="l")
plot(curve.to.prob(mult_rsv_curves(t1,t2,t3,t4,time.seq,beta0,beta1)),type="l")

##Likelihood function
with(OC,{
  dbinom(x=cases,size=num_infant,prob=curve.to.prob(mult_rsv_curves(t1,t2,t3,t4,time.seq,beta0,beta1)),
         log=TRUE)
})
