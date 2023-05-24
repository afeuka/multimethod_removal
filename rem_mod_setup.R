### Title: Function for calling various removal models
### Author: Abbey Feuka
### Date: 12MAY22
### Notes: specifies multi-method removal model using pre-determined data structure
# allows for one time-varying covariate and one site-varying covariate

#################################
rem.mod.setup <- function(dat, #data file, rem_sim.R generates example
                          seed, #random seed for initial values
                          timeCov=F, #time-varying covariate
                          spatCov=F, #spatial covariate
                          interact=F, #interaction between two covariates
                          REperiod=F,#random effect for period
                          REsite=F, #random effect for site
                          availCorrect=F, #correct for limited availability
                          np,#number of detection parameters (1 for single, nTypes for one per removal type)
                          nChains, #number of mcmc chains
                          detectFun, #detection function ()
                          abundFun,#abundance function ("pois" for Poisson and "nb" for negative binomial)
                          postPred=F,#generate posterior predictions
                          logProb=F, #calcualte log probabilities
                          forcast=F){ #forecast one period beyond data
  require(nimble)
  require(boot)
  
  #no valid model called
  nimbleMod = NA
  const = NA
  inits = NA
  
  #read in data
  suppressMessages(
    source("C:/Users/Abigail.Feuka/OneDrive - USDA/bTB/Spatial_Risk_Assessment/deerabundance/nimble/Functions/Removal data prep.R")
  )
  #indexing 
  M <- dat$nSites
  J <- dat$nRemovals
  time <- dat$nPeriods
  
  set.seed(seed)
  #initial values
  N <- lambda <- y.pred <- mu <- theta <- list()
  for(i in 1:nChains){
    if(forcast==T){
      lambda[[i]] <- matrix(log(rnorm((time+1)*M,3,0.1)),time+1,M)
      mu[[i]] <- matrix(exp(rnorm((time+1)*M,1,0.1)),time+1,M)
    } else {
      lambda[[i]] <- matrix(log(rnorm((time)*M,3,0.1)),time,M)
      mu[[i]] <- matrix(exp(rnorm((time)*M,1,0.1)),time,M)
    }
    # N[[i]] <- colSums(rem.dat$hca)+1000
    if(forcast==T){
      N[[i]] <- matrix(rep(rpois(1,10000),(time+1)*M),
                       nrow=time+1,ncol=M)
    } else {
      N[[i]] <- matrix(rep(rpois(1,10000),time*M),
                             nrow=time,ncol=M)
      }
    if(forcast==T){
      theta[[i]] <- array(runif(M*J*(time+1),0.001,0.2),dim=c(J,time+1,M))
      y.pred[[i]] <- array(rpois(M*J*(time+1),rpois(1,200)),dim=c(J,time+1,M))
    } else {
      theta[[i]] <- array(runif(M*J*(time),0.001,0.2),dim=c(J,time,M))
      y.pred[[i]] <- array(rpois(M*J*time,rpois(1,200)),dim=c(J,time,M))}
  }
  
  if(detectFun=="catcheff"){
    p <- matrix(rbeta(nChains*J,1,1),nChains,np)
  } else if(detectFun=="logit"){
    alpha <- cbind(rnorm(nChains,0.5,0.1),rnorm(nChains,-0.5,0.1))#matrix(rnorm(2*nChains,0,1),nChains,2)
  }
  
  beta <- matrix(rnorm(3*nChains,0.15,0.1),nChains,3)
  nu <- rnorm(nChains,0.15,0.1)
  taunu <- rgamma(nChains,2,0.1)
  eta <- rnorm(nChains,0.15,0.1)
  taueta <- rgamma(nChains,2,0.1)
  r <- rpois(nChains,30)
  mup <- rnorm(nChains,0,0.1)
  taulam <- exp(rnorm(nChains,-2,0.2))
  tautheta <- exp(rnorm(nChains,-2,0.2))
  taup <- exp(rnorm(nChains,-2,0.2))

  ### Errors 
  if(is.null(detectFun)){
    print("Select a detection function.")
  } 
  ### Models
  nimbleMod <- eval(str2lang(paste("nimbleCode({",
                                   if(detectFun=="catcheff"){
                                     if(np>1){
                                       "
          mup ~ dnorm(0,1)
          taup ~ dgamma(0.1,0.1)
          for(j in 1:np){
          logit(p[j]) ~ dnorm(mup,taup)}
        "
                                     } else if(np==1){
                                       "phi ~ dbeta(1,1)"
                                     }
                                   } else if(detectFun=="logit"){
                                     "
        for(i in 1:2){
        alpha[i] ~ dlogis(0,1)}"},
                                   if((REperiod==T | REsite==T) & timeCov==F & spatCov==F){ #random effects, no covs
                                     "
        beta ~ dlogis(0,1)"
                                   } else if(REperiod==F & REsite==F){ #no random effects, covs
                                     if(timeCov==T & spatCov==F){
                                       "
        for(i in 1:2){
          beta[i] ~ dlogis(0,1)}"
                                     } else if(timeCov==F & spatCov==T){
                                       "
        for(i in 1:2){
          beta[i] ~ dlogis(0,1)}"
                                     } else if (timeCov==T & spatCov==T){
                                       "
        for(i in 1:3){
        beta[i] ~ dlogis(0,1)}"
                                     }
                                   }else if((REperiod==T | REsite==T) & (timeCov==T | spatCov==T)){ #random effects, covs
                                     if(timeCov==T & spatCov==F){
                                       "
        for(i in 1:2){
          beta[i] ~ dlogis(0,1)}"
                                     } else if(timeCov==F & spatCov==T){
                                       "
        for(i in 1:2){
          beta[i] ~ dlogis(0,1)}"
                                     } else if (timeCov==T & spatCov==T){
                                       "
        for(i in 1:3){
        beta[i] ~ dlogis(0,1)}"
                                     }
                                   },
        "
        r ~ dunif(0,100)         
        taulam ~ dgamma(0.1,0.1)"
        ,
                                   if(REsite==T & timeCov==F & spatCov==F){
                                     "
        taunu ~ dgamma(0.1,0.1)
        
        for(i in 1:M){
          nu[i] ~ dnorm(beta,taunu)
          nu.star[i] <- nu[i] - beta}"
                                   } else if(REsite==T) {
                                     "
        taunu ~ dgamma(0.1,0.1)
        
        for(i in 1:M){
          nu[i] ~ dnorm(beta[1],taunu)
          nu.star[i] <- nu[i] - beta[1]}"
                                   },
                                   if(REperiod==T & timeCov==F & spatCov==F){
                                       "
                                       taueta ~ dgamma(0.1,0.1)
                                       for(t in 1:time){
                                         eta[t] ~ dnorm(beta,taueta)
                                         eta.star[t] <- eta[t] - beta}
                                       "
                                   } else if(REperiod==T){
                                       "
                                       taueta ~ dgamma(0.1,0.1)
                                       for(t in 1:time){
                                         eta[t] ~ dnorm(beta[1],taueta)
                                         eta.star[t] <- eta[t] - beta[1]}
                                       "
                                   },
                                   "
       for(i in 1:M){
          N[1,i] ~ dpois(10000)
          ", if(forcast==T){
            "
            for(t in 1:time){
            "} else if(forcast==F){
              "
            for(t in 1:(time-1)){
            "
            },
                                   if(REperiod==F & REsite==F & timeCov==F & spatCov==F){
                                       "
                                       mu[t,i] ~ dnorm(0,sd=0.5)"
                                   } else if(REperiod==F & REsite==F){ #no random effects
                                     if(timeCov==T & spatCov==T){
                                       "
          mu[t,i] <- beta[1] + beta[2]*sc[i] + beta[3]*tc[t]"
                                     } else if (timeCov==T & spatCov==F){
                                       "
          mu[t,i] <- beta[1] + beta[2]*tc[t]"
                                     } else if (timeCov==F & spatCov==T){
                                       "
          mu[t,i] <- beta[1] + beta[2]*sc[i]"
                                     }        
                                   } else if(REperiod==T & REsite==F){ #random effects year
                                     if(timeCov==T & spatCov==T){
                                       "
          mu[t,i] <- eta[t] + beta[2]*sc[i] + beta[3]*tc[t]"
                                     } else if (timeCov==T & spatCov==F){
                                       "
          mu[t,i] <- eta[t] + beta[2]*tc[t]"
                                     } else if (timeCov==F & spatCov==T){
                                       "
          mu[t,i] <- eta[t] + beta[2]*sc[i]"
                                     } else if (timeCov==F & spatCov==F){
                                       "
          mu[t,i] <- eta[t]"
                                     }
                                   } else if(REperiod==F & REsite==T){ #random effects site
                                     if(timeCov==T & spatCov==T){
                                       "
          mu[t,i] <- nu[i] + beta[2]*sc[i] + beta[3]*tc[t]"
                                     } else if (timeCov==T & spatCov==F){
                                       "
          mu[t,i] <- nu[i] + beta[2]*tc[t]"
                                     } else if (timeCov==F & spatCov==T){
                                       "
          mu[t,i] <- nu[i] + beta[2]*sc[i]"
                                     } else if (timeCov==F & spatCov==F){
                                       "
          mu[t,i] <- nu[i]"
                                     }
                                   } else if(REperiod==T & REsite==T){ #random effects year and site
                                     if(timeCov==T & spatCov==T){
                                       "
          mu[t,i] <- nu[i] + eta[t] + beta[2]*sc[i] + beta[3]*tc[t]"
                                     } else if (timeCov==T & spatCov==F){
                                       "
          mu[t,i] <- nu[i] + eta[t] + beta[2]*tc[t]"
                                     } else if (timeCov==F & spatCov==T){
                                       "
          mu[t,i] <- nu[i] + eta[t] + beta[2]*sc[i]"
                                     } else if (timeCov==F & spatCov==F){
                                       "
          mu[t,i] <- nu[i] + eta[t]"
                                     }
                                   },
                                   if(abundFun=="pois"){
                                     "
            log(lambda[t,i]) ~ dnorm(mu[t,i],taulam)                        
            delta[t,i] <- lambda[t,i]*(N[t,i]-sum(y[1:J,t,i]))
            N[t+1,i] ~ dpois(delta[t,i])
                                     }"
                                   } else if(abundFun=="nb"){
                                     "
        log(lambda[t,i]) ~ dnorm(mu[t,i],taulam) 
        delta[t,i] <- lambda[t,i]*(N[t,i]-sum(y[1:J,t,i]))
        nb.p[t,i] <- r/(delta[t,i]+r)
        N[t+1,i] ~ dnegbin(prob=nb.p[t,i],size=r)
            }"
        },
        if(forcast==T){
          "
        for(t in 1:(time+1)){
          "} else if(forcast==F){
            "
        for(t in 1:time){
          "},
        "
          for(j in 1:J){
            ",
                                   if(detectFun=="catcheff"){
                                     if(np==J){
                                       "
          theta[j,t,i] <- 1-pow((1-p[j]),e[j,t,i])"
                                     } else if(np==1){
                                       "
          theta[j,t,i] <- 1-pow((1-p),e[j,t,i])"}                             
                                   } else if(detectFun=="logit"){
                                     "
        logit(theta[j,t,i]) ~ alpha[1] + alpha[2]*e[j,t,i]"},
                                   if(availCorrect==TRUE){
                                     "
        pip[j,t,i] <- avail_fun(removal=j,
                                gamma=gamma[j,i],
                                gamma.past=gamma[1:(j-1),i],
                                theta=theta[j,t,i],
                                theta.past=theta[1:(j-1),t,i])"
                                   } else {
                                     "
        pip[j,t,i] <- theta[j,t,i]"},
         if(postPred==TRUE){
          "
           y.pred[j,t,i] ~ dbin(prob=pip[j,t,i],size=N[t,i])"
         },
        "
        }", #removals
        "}",#years
          "
          for(t in 1:time){
            for(j in 1:J){
              y[j,t,i] ~ dbin(prob=pip[j,t,i],size=N[t,i])
            }
          }",
        "}",#sites
                                   if(logProb==T){
                                     "
      sumLogProbY ~ dnorm(0,1) 
      sumLogProbY.pred ~ dnorm(0,1)"
                                   },
                                   "})"))#end of nimbleCode and str2lang
  ) #end of eval
  
  # assign constants 
  const <- list(J = dat$nRemovals,
                time = dat$nYears,
                M = dat$nSites,
                e = dat$eff,
                A = dat$A,
                if(detectFun=="catcheff") np = dim(dat$hca)[1]
  )
  
  const <- const[!sapply(const,is.null)]
  
  if(timeCov==T | spatCov==T){
    names(const) <- c("J","time","M","e","A")
  }
  if(timeCov==F & spatCov==F){
    names(const) <- c("J","time","M","e","A")
  }
  if(detectFun=="catcheff"){
    names(const)[length(names(const))] <- "np"
  }
  
  init.names <- c("N",
                  if(detectFun=="catcheff"){
                    "p"
                  } else if(detectFun=="logit"){
                    "alpha"
                  },
                  if(detectFun=="catcheff"){
                    "theta"
                  },
                  if(detectFun=="catcheff"){
                    "tautheta"
                  },
                    "lambda"
                  ,
                    "mu"
                  ,
                    "taulam"
                  ,
                  if(timeCov==T & spatCov==T){
                    "beta"
                  } else if(timeCov==T & spatCov==F){
                    "beta"
                  } else if(timeCov==F & spatCov==T){
                    "beta"
                  } else if(REsite==T | REperiod==T){
                    "beta"
                  },
                  if(REperiod==T){
                    "eta"},
                  if(REperiod==T){
                    "taueta"},
                  if(REsite==T){
                    "nu"},
                  if(REsite==T){
                    "taunu"},
                  if(abundFun=="nb"){
                    "r"
                  },
                  if(postPred==T){
                    "y.pred"
                  },
                  if(logProb==T){
                    "sumLogProbY"},
                  if(logProb==T){
                    "sumLogProbY.pred"
                  })    

  inits <- list()
  for(i in 1:nChains){
    inits[[i]] <- list(
      N[[i]],
      if(detectFun=="catcheff"){
        p[i,]
      } else if(detectFun=="logit"){
        alpha[i,]
      },
      if(detectFun=="catcheff"){
        theta[[i]]
      },
      if(detectFun=="catcheff"){
        tautheta[i]
      },
        lambda[[i]]
      ,
        mu[[i]]
      ,
        taulam[i]
      ,
      if(timeCov==T & spatCov==T & interact==T){
        beta[i,]
      } else if(timeCov==T & spatCov==T & interact==F){
        beta[i,1:3]
      } else if(timeCov==T & spatCov==F){
        beta[i,1:2]
      } else if(timeCov==F & spatCov==T){
        beta[i,1:2]
      } else if(REsite==T | REperiod==T){
        beta[i,1]
      },
      if(REperiod==T){
        rep(eta[i],time)},
      if(REperiod==T){
        taueta[i]},
      if(REsite==T){
        rep(nu[i],M)},
      if(REsite==T){
        taunu[i]},
      if(abundFun=="nb"){
        r[i]},
      if(postPred==T){
        y.pred[[i]]
      },
      if(logProb==T){
        0
      },
      if(logProb==T){
        0
      }
    )
    inits[[i]] <- inits[[i]][!sapply(inits[[i]],is.null)]
    names(inits[[i]]) <- init.names
  }
  
  
  list(nimbleMod = nimbleMod, const = const, inits = inits)
  
} #end of rem.mod.cond()

