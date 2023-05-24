### Title: Fitting removal models to estimate deer abundance
### Author: Abbey Feuka
### Date: 25MAY22
### Notes: 
############################
library(nimble)
library(MCMCvis)

#simulated data
seed<-1
nperiods <- 10
nsites <- 30
nremovals <- 5
gam <- matrix(rbeta(nremovals*nsites),nremovals,nsites)
eff <- array(scale(rpois(nremovals*nsites*nperiods,10)),
             dim=c(nremovals,nperiods,nsites))
area <- rpois(nsites,30)
tc <- rnorm(nperiods,0,1)
sc <- rnorm(nsites,0,1)
source("./Functions/rem_sim.R")
sim <- rm.sim(seed=seed,tc=tc,sc=sc,gam=gam,eff=eff,
              A=area,abund.dist="pois",det.fun="logit",
              capPerUnitEff=c(-4,1),beta=c(-0.5,0.4,-0.2),
              r=30,REperiod=F,REsite=F,
              avail.correct=T,plots=T)

modsetDat <- list(y=sim$y,
                  tc=tc,
                  sc=sc,
                  gamma=gam)

# setup model structure fit to all data structures
niter <- 200000
burnProp <- 0.5
n.chains <- 3
seed <- 1:n.chains

source("./Functions/rem_mod_setup.R") # function for specifying model structure
modset <- rem.mod.setup(dat=sim, seed=seed,
                        timeCov=T, spatCov=F, REperiod=F, REsite=F,
                        detectFun="catcheff", abundFun = "pois",
                        availCorrect=T, np=nremovals,
                        nChains=n.chains, postPred =T,
                        logProb=T,forecast=F)

monitors <- c(names(modset$inits[[1]]),"pip")
modset$nimbleMod

#####################
### Fit model - single chain
source("./Functions/run.mcmc.allcode.R") # function for fitting model
mcmc.out <- run_MCMC_allcode(seed=seed,
  mod=modset$nimbleMod,
  dat=modsetDat,
  const=modset$const,
  inits=modset$inits[[1]],
  niter=niter,standard=T,
  params=sapply(1:6,function(i)paste0("N[1, ",i,"]")),
  sampler.typ="RW",monitors=monitors,nchains=n.chains,
  burnProp=burnProp,compile=T,run=T,log.probs=T,
  waic=F,control=list(scale=10))

samples <- mcmc.out$mcmc.out
samples <- rbind.data.frame(samples$chain1,samples$chain2,samples$chain3)
samples$chain <- sort(rep(1:n.chains,nrow(samples)/n.chains))

#multiple chains
# source("./Functions/fit.rem.mod.R")
# mcmc.out <- fit.rem.mod(seed=seed,modsetDat = modsetDat,modset=modset,
#                         monitors=monitors,n.chains=n.chains,niter=niter,
#                         burnProp = burnProp,standard=T)

# mod <- mcmc.out$mod
# cmod <- mcmc.out$Cmod
# mcmc <- mcmc.out$mcmc.uncomp
# cmcmc <- mcmc.out$Cmcmc
# samples<-list()
# for(k in 1:n.chains){
#   samples[[k]] <- MCMCchains(mcmc.out$parSamples[[k]]$mcmc.out$samples,
#                              params=c("beta","alpha[1]","alpha[2]",
#                                       sapply(1:11,function(t){sapply(1:6,function(i){paste0("N[",t,", ",i,"]")})}),
#                                       sapply(1:10,function(t){sapply(1:6,function(i){paste0("lambda[",t,", ",i,"]")})})),
#                              #sapply(1:11,function(t){sapply(1:6,function(i){paste0("pip[",t,", ",i,"]")})})),
#                              ISB=F)
# }

## traceplots
source("./Functions/traceplot_fun.R")
samples <- rbind.data.frame(samples$chain1,samples$chain2,samples$chain3)
samples$chain <- sort(rep(1:3,niter*burnProp))
samples$idx <- rep(1:(niter*burnProp),3)
samples <- samples %>%
  pivot_longer(cols=1:(ncol(samples)-2),names_to = "parameter",values_to="value")

# bayesian p-values using deviance
dev.pred <- -2*samples[,'sumLogProbY.pred']
dev.y <- -2*samples[,'sumLogProbY']
ind <- dev.pred > dev.y
mean(ind)

par(mfrow=c(1,1))
matplot(MCMCchains(samples,
                   params = "beta",ISB=T),
        typ="l",ylab="betas")
abline(h=sim$beta,col=1:length(sim$beta))

par(mfrow=c(1,1))
matplot(MCMCchains(samples,params="alpha"),typ="l")
abline(h=sim$alpha,col=1:length(sim$alpha))

hist(MCMCchains(samples,params="pip"),freq=F)
hist(sim$pip,col=alpha("red",0.5),add=T,freq=F)

cor(cbind(samples[,"beta[1]"],samples[,"beta[2]"],
          samples[,"alpha[1]"],samples[,"alpha[2]"]))

## posterior prediction
source("./Functions/postpred_checks.R")
pp <- post.pred(y.pred=samples[,grep("y.pred",colnames(samples))],rem.dat=rem.dat)
