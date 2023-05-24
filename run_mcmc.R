### Title: Fit removal nimbleModel to data
### Author: Abbey Feuka
### Date: 13MAY22
### Notes: Useful for parallelization, outputs mcmc samples and other mcmc objects
##########################

run_MCMC_allcode <- function(seed, #random seed for mcmc
                             mod, #nimbleModel generated from rem_mod_setup
                             dat, #datafile from rem_sim or real data
                             const, #list of constants from rem_mod_setup
                             inits, #list of initial values from rem_mod_setup
                             niter, #number of mcmc iterations
                             standard=TRUE, #standard MCMC algorithm determined by nimble
                             params=NA, #parameters to block, if desired
                             sampler.typ=NA, #specify if blocking groups of parameters: "AF_slice" or "RW_block"
                             control=NA, #control parameters for MCMC, if desired
                             monitors, #parameters to monitor
                             burnProp, #proportion of chains as burn in
                             nchains=1, #number of chains
                             compile=FALSE, #compile nimble model and MCMC
                             run=FALSE, #runs mcmc algorithm/fits model
                             log.probs=FALSE, #calculate log probabilities of data
                             Cmcmc=NA, #insert pre-made compiled mcmc if have one
                             waic=F) { #calculate waic
  require(nimble)

  source("./Functions/avail_fun.R")
  source("./Functions/storing_log_probs.R")

  if(is.na(Cmcmc)){
    mod <- nimbleModel(code = mod,
                       data = dat,
                       constants = const,
                       inits = inits)
    
    mcmc.conf <- configureMCMC(mod, enableWAIC=waic, print = TRUE) #default MCMC configuration
    # mcmc.conf <- configureMCMC(mod, enableWAIC=waic, print = TRUE, 
    #                            autoBlock=T) #tries different parameter blockings (RW_block) to increase efficency
    mcmc.conf$addMonitors(monitors)
    #storing log probs for data and posterior predictions
    if(log.probs==TRUE){
      configureStoreLogProb(mcmc.conf, mod, 'sumLogProbY',
                            mod$getNodeNames(dataOnly = TRUE))
      configureStoreLogProb(mcmc.conf, mod, 'sumLogProbY.pred',
                            mod$getNodeNames()[grep("y.pred",mod$getNodeNames())])
    }
    mcmc.uncomp <- buildMCMC(mcmc.conf) #uncompiled MCMC
  } else {
    mcmc.conf <- NA
    mcmc.uncomp <- NA
  }
  
  if(compile==TRUE){
    Cmod <- compileNimble(mod) #compiled model
  } else {
    Cmod <- NA
  }

  if(standard==TRUE){
    if(compile==TRUE){
      Cmcmc <- compileNimble(mcmc.uncomp, project=mod) #compiled mcmc
    } 
  } else if(standard==FALSE){
    # remove samplers
      if(!(sampler.typ %in% c("RW_block","AF_slice"))){
        for(i in 1:length(params)){
          mcmc.conf$removeSamplers(params[i])
          mcmc.conf$addSampler(target = params[i] ,type = sampler.typ)
        }
      } else if(sampler.typ %in% c("RW_block","AF_slice")){
        for(i in 1:length(params)){
          mcmc.conf$removeSamplers(params[i])
        }
        if(!is.na(control)){
          mcmc.conf$addSampler(target = params ,type = sampler.typ,
                               control = control)
        } else {
          mcmc.conf$addSampler(target = params ,type = sampler.typ)
        }

      }
    
    if(is.na(Cmcmc)){
      mcmc.samplers <- mcmc.conf$printSamplers()
      mcmc.conf <- buildMCMC(mcmc.conf)
      # need to reset the nimbleFunctions in order to add the new MCMC
      Cmcmc <- compileNimble(mcmc.conf, resetFunctions = TRUE)
    }
  }
  
  if(run==TRUE){
    mcmc.out <- runMCMC(Cmcmc, niter = niter, nburn = niter*burnProp,
                       setSeed = seed, nchains = nchains, WAIC = waic)
  } else {
    mcmc.out <- NA
  }

  out <- list(mcmc.out=mcmc.out,
              mod=mod,mcmc.uncomp=mcmc.uncomp,
              mcmc.conf=mcmc.conf,
              Cmod=Cmod,Cmcmc=Cmcmc)
  return(out)
}
