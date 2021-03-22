###########################
## Libraries that I need ##
###########################
library(tidyverse)
library(reshape)

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

rsv_curve <- function(t1,t2,t3,t4,time, beta0,theta){
  out <- rep(beta0,length(time))
  out[t1<time & time <=t2] <- beta0 + theta*(time[t1<time & time <=t2]-t1)/(t2-t1)
  out[t2<time & time <=t3] <- beta0 + theta
  out[t3<time & time <=t4] <- beta0 + theta - theta*(time[t3<time & time <=t4]-t3)/(t4-t3)
  return(out)
}
#qplot(1:24,rsv_curve(6,11,14,17,1:24,-2,1),geom="line")

mult_rsv_curves <- function(t1,t2,t3,t4,time,beta0,theta){
  n.years <- length(theta)
  return(c(sapply(1:n.years,function(x){rsv_curve(t1[x],t2[x],t3[x],t4[x],time,beta0,theta[x])}))) #beta0 constant over years?
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
# the.cycle <- 10
# OC %>% dplyr::filter(cycle==the.cycle) %>% ggplot(aes(x=1:24,y=cases/num_infant))+geom_line()
# OC <- OC %>% dplyr::filter(cycle==the.cycle)
num.years <- length(unique(OC$cycle))

# #######################################
# ## Read in the Data from Salt Lake ##
# #######################################
# which(state=='utah' & county=='utah')
# 
# SL <- RSV[[2769]]$rsv.df # SLC
# SL <- SL[which(SL$Strt=="2003-07-01"):which(SL$End == "2013-06-30"),-1]
# SL <- within(SL,
#             {cycle <- rep(1:(nrow(SL)/24), each=24)
#             })
# num.years <- max(SL$cycle)


###################################################################
## Choose starting values for xi1[1:num.years],xi2[1:num.years], ##
## xi3[1:num.years],xi4[1:num.years],beta0,beta1[1:num.years]    ##
###################################################################
xi <- matrix(-1,nrow=num.years,ncol=4)
delta <- t(apply(xi,1,xi.to.delta))
chg.pts <- t(apply(delta,1,delta.to.t))
beta0 <- -10
theta <- (OC %>% dplyr::group_by(cycle) %>% summarize(theta=log(max(cases/num_infant)/(1-max(cases/num_infant)))-beta0))[['theta']]
current.probs <- curve.to.prob(mult_rsv_curves(chg.pts[,1],chg.pts[,2],chg.pts[,3],chg.pts[,4],1:24,beta0,theta))
parameters <- c(as.vector(xi), beta0, log(theta))
# plot(x=1:length(current.probs),y=current.probs,type="l",ylim=c(0,.03))
# points(OC$cases/OC$num_infant)

#########################
## MCMC specifications ##
#########################
source("./AMCMCUpdate.R")
n.parameters <- length(xi)+length(beta0)+length(theta)
amcmc <- list(mn=matrix(0,nrow=n.parameters,ncol=1),var=matrix(0,nrow=n.parameters,ncol=n.parameters))
amcmc.it <- 100
n.it <- 1000
burn <- 10000
thin <- 50
tot.it <- burn+thin*n.it
kp.seq <- seq(burn+thin,tot.it,by=thin)
kp <- 0

############
## Priors ##
############
pri.mn <- c(rep(0,length(xi)),-10,5)
pri.sd <- rep(1,n.parameters)

#################################
## Matrices to hold MCMC draws ##
#################################
chgpt.draws <- array(0,dim=c(num.years,4,n.it))
beta0.draws <- rep(0,n.it)
theta.draws <- matrix(0,nrow=n.it,ncol=num.years)

##################
## Run the MCMC ##
##################
for(it in 1:tot.it){
  
  ## Propose new parameters for all years
  if(it<=amcmc.it){
    prop.var <- (0.001^2)*diag(n.parameters)
  } else {
    prop.var <- (2.4^2/n.parameters)*((.001^2)*diag(n.parameters)+ amcmc$var)
  }
  prop.parameters <- c(xi,beta0,log(theta))+t(chol(prop.var))%*%rnorm(n.parameters)
  #prop.parameters[length(xi)+1] <- beta0
  #prop.parameters[length(xi)+1+(1:num.years)] <- log(theta)
  prop.xi <- matrix(prop.parameters[1:length(xi)],nrow=num.years,ncol=4)
  prop.delta <- t(apply(prop.xi,1,xi.to.delta))
  prop.chgpt <- t(apply(prop.delta,1,delta.to.t))
  prop.beta0 <- prop.parameters[length(xi)+1]
  prop.theta <- exp(prop.parameters[length(xi)+1+(1:num.years)])
  
  ## Calculate Proposed probabilities
  prop.probs <- curve.to.prob(mult_rsv_curves(prop.chgpt[,1],prop.chgpt[,2],prop.chgpt[,3],prop.chgpt[,4],1:24,prop.beta0,prop.theta))
  
  ## Check Metropolis Hastings acceptance probability c(as.vector(xi), beta0, beta1)
  prior.star <- dnorm(prop.parameters, mean=pri.mn, sd=pri.sd, log=TRUE)
  prior <- dnorm(parameters, mean=pri.mn, sd=pri.sd, log=TRUE)
  
  numerator <- sum(dbinom(OC$cases,OC$num_infant,prob=prop.probs,log=TRUE)) + sum(prior.star)
  denominator <- sum(dbinom(OC$cases,OC$num_infant,prob=current.probs,log=TRUE)) + sum(prior)
  r <- numerator - denominator
  #print(r)
  
  ## Keep or throw proposed values
  if(log(runif(1))<r){
    parameters <- prop.parameters
    xi <- prop.xi
    delta <- prop.delta
    chg.pts <- prop.chgpt
    beta0 <- prop.beta0
    theta <- prop.theta
    current.probs <- prop.probs
  }
  
  ## Update the AMCMC proposal variance
  amcmc <- AMCMC.update(parameters,amcmc$mn,amcmc$var,it) 
  
  ## Keep the draw if necessary
  if(it%in%kp.seq){
    kp <- kp+1
    chgpt.draws[,,kp] <- chg.pts
    beta0.draws[kp] <- beta0
    theta.draws[kp,]<- theta
  }
  
  
} #end M
#plot(theta.draws[,2],type="l")
chgpt.est <- apply(chgpt.draws, c(1,2), mean)
theta.est <- apply(theta.draws, 2, mean)
beta0.est <- mean(beta0.draws)

the.years <- 3:7
plot.prob <- curve.to.prob(mult_rsv_curves(chgpt.est[the.years,1], 
                                           chgpt.est[the.years,2], 
                                           chgpt.est[the.years,3], 
                                           chgpt.est[the.years,4],
                                           1:24, 
                                           beta0.est,
                                           theta.est[the.years]))
obs.probs <- OC$cases[OC$cycle%in%the.years]/OC$num_infant[OC$cycle%in%the.years]
plot(x=1:length(plot.prob),y=plot.prob,type="l",ylim=range(obs.probs,plot.prob))
points(obs.probs,pch=19)
