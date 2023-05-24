### Title: Simulating removal data with real covariates
### Author: Abbey Feuka
### Date: 05MAY22
### Notes: 
##########################
rm.sim <- function(seed=1,#random seed for simulation
                   tc, #time-varying covariate
                   sc, #spatially varying covariate
                   gam, #matrix of proportion available where nrow=number of removals and ncol=site
                   eff, #array of effort with dim=c(n survey types, nPeriods, nsites)
                   A, #area of sites
                   abund.dist, #abundance distribution, poisson "pois" or negative binomial "nb"
                   det.fun, #choose detection function ("logit" for logit-linear regression and "catcheff" for catch-effort curve)
                   capPerUnitEff, #capture rate per unit effort (vector minimum of 2 for logit specification and scalar for catcheff)
                   beta=c(1,1,1), #coefficients on growth rate c(intercept,spatial covariate,time covariate)
                   r=NA, #size parameter for negative binomial distribution, if using
                   REperiod=FALSE,REsite=FALSE, #random effects for site or period
                   avail.correct=FALSE, #correction for incomplete availability
                   plots=FALSE){ #produce plots of simulated abundance
  set.seed(seed)
  
  source("./Functions/avail.fun.R")
  #errors
  if(sum(is.na(capPerUnitEff))!=0 & det.fun=="logit"){
    print("Assign capPerUnitEff values")
  }
  if(sum(is.na(phi))!=0 & det.fun=="catcheff"){
    print("Assign capPerUnitEff value")
  }

  #index parameters 
  nRemovals <- dim(gam)[1]
  nPeriods <- length(tc) #years
  nSites <- length(sc) #sites
  
  #density simulation
  N <- y <- pip <- theta <- array(rep(matrix(NA,nRemovals,nPeriods),nSites),dim=c(nRemovals,nPeriods,nSites))
  lambda <- matrix(NA,nrow=nPeriods,ncol=nSites)
  if(abund.dist=="pois"){
    N <- lambda <- matrix(NA,nrow=nPeriods,ncol=nSites)
  }

  eta <- rnorm(nPeriods,0.2,0.1)
  nu <- rnorm(nSites,0.3,0.1)
  for(i in 1:nSites){
    for(t in 1:nPeriods){
      if(REperiod==F & REsite==F){
        lambda[t,i] <- exp(beta[1] + beta[2]*sc[i] + beta[3]*tc[t])
      } else if(REperiod==F & REsite==T){
        lambda[t,i] <- exp(beta[1] + beta[2]*sc[i] + beta[3]*tc[t] + nu[i])
      } else if(REperiod==T & REsite==F){
        lambda[t,i] <- exp(beta[1] + beta[2]*sc[i] + beta[3]*tc[t] + eta[t])
      } else if(REperiod==T & REsite==T){
        lambda[t,i] <- exp(beta[1] + beta[2]*sc[i] + beta[3]*tc[t] + 
                            nu[i] + eta[t])
      }
    }
    
    for(t in 1:nPeriods){
      if(abund.dist=="pois"){
          if(t==1){
            N[t,i] <- rpois(1,10000)
          } else {
            delta[t-1,i] <- lambda[t-1,i]*(N[t-1,i]-sum(y[1:nRemovals,t-1,i]))
            N[t,i] <- rpois(1,delta[t-1,i])
          }
      } else if(abund.dist=="nb"){
          mu <- lambda[t,i]*A[i]
          nb.p <- r/(mu+r)
          N[1,t,i] <- rnbinom(1,prob=nb.p,size=r)
      }
      for(j in 1:nRemovals){
        if(det.fun=="logit"){
          theta[j,t,i] <- inv.logit(capPerUnitEff[1] + capPerUnitEff[2]*eff[j,t,i])
        } else if(det.fun=="catcheff"){
          theta[j,t,i] <- 1-(1-capPerUnitEff[j])^eff[j,t,i]
        }
        if(avail.correct==T){
          pip[j,t,i] <- avail_fun(removal=j,
                                  gamma=gam[j,i],
                                  gamma.past=gam[1:(j-1),i],
                                  theta=theta[j,t,i],
                                  theta.past=theta[1:(j-1),t,i])
        } else {
          pip[j,t,i] <- theta[j,t,i]
        }
        if(abund.dist=="nb"){
          if(j==1){
            y[1,t,i] <- rbinom(1,N[1,t,i],pip[1,t,i])
          } else if(j>1){
            N[j,t,i] <- N[j-1,t,i] - y[j-1,t,i]
            y[j,t,i] <- rbinom(1,N[j,t,i],pip[j,t,i])
          }
        } else if(abund.dist=="pois"){
          y[j,t,i] <- rbinom(1,N[t,i],pip[j,t,i])
        }
      } #nRemovals loop
      } #nPeriods loop
      }#nSites loop
  
  
  if(plots==TRUE){
    if(abund.dist=="nb"){
      N.plot <- matplot(t(N[,,1]),typ="l",ylab="Alcona DMU Abundance")
      legend("topleft",legend=sapply(1:nRemovals,function(i){paste("Hunt",i)}),
             col=1:nRemovals,lty=1:nRemovals,horiz = T,xpd="o",bty="n")
    } else if(abund.dist=="pois"){
      N.plot <- matplot(N,typ="l",ylab="Alcona DMU Abundance")
      legend("topleft",legend=sapply(1:nRemovals,function(i){paste("Hunt",i)}),
             col=1:nRemovals,lty=1:nRemovals,horiz = T,xpd="o",bty="n")
    }
    y.plot <- matplot(t(y[,,1]),typ="l",ylab="Alcona Harvest")
    legend(x=c(0,10),y=c(-50,-100),legend=sapply(1:nRemovals,function(i){paste("Hunt",i)}),
           col=1:nRemovals,lty=1:nRemovals,horiz = T,xpd="o",bty="n")
    hist(pip)
  } else {
    N.plot <- NA
    y.plot <- NA
  }
  
  return(list(y=y,N=N,A=A,pip=pip,beta=beta,r=r,capPerUnitEff=capPerUnitEff,
              eta=eta,nu=nu,
              lambda=lambda,gam=gam,eff=eff,nRemovals=nRemovals,
              nSites=nSites,nRemovals=nRemovals,nPeriods=nPeriods,
              y.plot=y.plot,N.plot=N.plot))
  
}

