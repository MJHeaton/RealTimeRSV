###########################
## Libraries that I need ##
###########################
require(ggplot2)
library(reshape)
library(magrittr) #for piping %>%

###############################
## Functions for Convenience ##
###############################
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
  n.years <- length(beta1)
  return(c(sapply(1:n.years,function(x){rsv_curve(t1[x],t2[x],t3[x],t4[x],time,beta0,beta1[x])}))) #beta0 constant over years?
}

curve.to.prob <- function(curve.val){
  return(plogis(curve.val))
}

##############################
## Read in the Data from OC ##
##############################
OC <- read.csv(file="theOC.csv")
OC <- OC[which(OC$Strt=="2003-07-01"):which(OC$End=="2013-06-30"),-1]
OC <- within(OC,{
  cycle <- rep(1:(nrow(OC)/24),each=24)
})
num.years <- max(OC$cycle)

###################################################################
## Choose starting values for xi1[1:num.years],xi2[1:num.years], ##
## xi3[1:num.years],xi4[1:num.years],beta0,beta1[1:num.years]    ##
###################################################################
xi <- matrix(-1,nrow=num.years,ncol=4)
delta <- t(apply(xi,1,xi.to.delta))
chg.pts <- t(apply(delta,1,delta.to.t))
beta0 <- -3
beta1 <- rep(.5,num.years)
current.probs <- curve.to.prob(mult_rsv_curves(chg.pts[,1],chg.pts[,2],chg.pts[,3],chg.pts[,4],1:24,beta0,beta1))
qplot(x=1:length(current.probs),y=current.probs,geom="line")

#########################
## MCMC specifications ##
#########################
source("~/Documents/Professional/Computation/RFunctions/AMCMCUpdate.R")
n.parameters <- prod(dim(xi))+length(beta0)+length(beta1)
amcmc <- list(mn=matrix(0,nrow=n.parameters,ncol=1),var=matrix(0,nrow=n.parameters,ncol=n.parameters))
amcmc.it <- 100
n.it <- 1000
burn <- 100
thin <- 1
tot.it <- burn+thin*n.it
kp.seq <- seq(burn+thin,tot.it,by=thin)

#################################
## Matrices to MCMC hold draws ##
#################################
chgpt.draws <- array(0,dim=c(num.years,4,n.it))
beta0.draws <- rep(0,n.it)
beta1.draws <- matrix(0,nrow=n.it,ncol=num.years)

##################
## Run the MCMC ##
##################
for(i in 1:tot.it){
  
  ## Propose new parameters for all years
  if(it<=amcmc.it){
    prop.var <- (0.001^2)*diag(n.parameters)
  } else {
    prop.var <- (2.4^2/n.parameters)*((.001^2)*diag(num.years)+params.amcmc$var)
  }
  prop.parameters <- c(xi,beta0,beta1)+t(chol(prop.var))%*%rnorm(n.parameters)
  prop.xi <- matrix(prop.parameters[1:length(xi)],nrow=num.years,ncol=4)
  prop.delta <- t(apply(prop.xi,1,xi.to.delta))
  prop.chgpt <- t(apply(prop.delta,1,delta.to.t))
  prop.beta0 <- prop.parameters[length(xi)+1]
  prop.beta1 <- prop.parameters[length(xi)+1+(1:num.years)]
  
  ## Calculate Proposed probabilities
  prop.probs <- curve.to.prob(mult_rsv_curves(prop.chgpt[,1],prop.chgpt[,2],prop.chgpt[,3],prop.chgpt[,4],1:24,prop.beta0,prop.beta1))
  
  ## Check Metropolis Hastings acceptance probability
  numerator <- sum(dbinom(OC$cases,OC$num_infant,prob=prop.probs,log=TRUE))+
    prior
  denominator <- fill.in.on.your.own
  r <- fill.in.on.your.own
  
  ## Update the AMCMC proposal variance
  
  ## Keep the draw if necessary
  if(it%in%kp.seq){
    kp <- kp+1
    chgpt.draws[,,kp] <- chg.pts
  }
  
  
} #end MCMC for-loop




