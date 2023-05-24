### Title: Function for plotting trace plots and cumulative mean plots
### Author: Abbey Feuka
### Date: 20SEP22
### Notes:

tracePlot <- function(samples, #mcmc samples generated from run_mcmc
                      param, #parameters to plot
                      nChains, #number of chains
                      cumPlot = FALSE){ #cumulative mean plot to assess convergence
  if(cumPlot==T){
    par(mfrow=c(1,2))
  } else {
    par(mfrow=c(1,1))
  }
  # if(is.list(samples)){
  #   #trace
  #   plot(samples[[1]][,param],typ="l",ylab=param)
  #   for(i in 2:nChains){
  #     points(samples[[i]][,param],typ="l",col=i)
  #   }
  #   if(cumPlot==T){
  #     #cumulative mean
  #     plot(cumsum(samples[[1]][,param])/
  #            (1:nrow(samples[[1]])),typ="l",ylab=paste0("mean(",param,")"))
  #     for(i in 2:nChains){
  #       points(cumsum(samples[[i]][,param])/
  #                (1:nrow(samples[[i]])),
  #              typ="l",col=i)
  #     }
  #   }
  # }else {
  if("chain"%in%colnames(samples)){
    #trace
    require(tidyverse)
    df <- samples %>% filter(parameter==param)
    df$chain <- as.factor(df$chain)
    ggplot(df,aes(y=value,x=idx,col=chain)) + geom_line() +
      theme(panel.background = element_blank())
    # if(cumPlot==T){
    #   #cumulative mean
    #   plot(cumsum(samples[,param])/
    #          (1:nrow(samples)),typ="l",ylab=paste0("mean(",param,")"))
    # }
  } else {
    #trace
    plot(samples[,param],typ="l",ylab=param)
    points(samples[,param],typ="l",col=i)
    if(cumPlot==T){
      #cumulative mean
      plot(cumsum(samples[,param])/
             (1:nrow(samples)),typ="l",ylab=paste0("mean(",param,")"))
    }
  }
  #}
}