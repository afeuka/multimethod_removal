### Title: Posterior checks on abundance model
### Author: Abbey Feuka
### Date: 19MAY22
### Notes: Uses outputs from run_mcmc/rem_mod_script to generate posterior predictive plots compared to data
############################
post.pred <- function(y.pred, #posterior predictions from run_mcmc
                      dat, #full data set
                      removal.names=NA, #names of removals, if named, in order appearing in dat
                      # mod,  #mcmc.out$mod
                      # mcmc, #mcmc.out$mcmc.conf
                      sites.sf=NA, #sf object of sites, linked to index numbers used in model
                      legend=F){ #plot legend
  require(nimble)
  require(tidyverse)

# if(identical(colnames(samples), 
#           mod$expandNodeNames(mcmc$mvSamples$getVarNames()))==FALSE){
#   ord <- numeric(ncol(samples))
#   for(i in seq_len(length(colnames(samples)))){
#     ord[i] <- which(mod$expandNodeNames(mcmc$mvSamples$getVarNames())[i]==colnames(samples))
#   }
#   samples <- samples[,ord] #reorder to match
# }

# identical(colnames(samples), 
#           mod$expandNodeNames(mcmc$mvSamples$getVarNames()))



## plotting pp
# pp <- samples[,mod$getNodeNames()[grep("y.pred",mod$getNodeNames())]]
pp <- y.pred
n.preds <- ncol(pp)

if(is.na(removal.names)){
  removal.names <- 1:dat$nRemovals
}
y <- dat$y
#dataframe indexing data/posterior preds vectors
ppidx <- substr(colnames(pp),8,20)
ppidx <- strsplit(ppidx,",")
for(i in 1:length(ppidx)){
  ppidx[[i]] <- as.data.frame(t(ppidx[[i]]))
}
ppidx <- do.call("rbind",ppidx)
for(i in 1:nrow(ppidx)){
  ppidx[i,3] <- as.numeric(substr(ppidx[i,3],2,2))
}
colnames(ppidx) <- c("season","year","dmu")

#point estimates and cred intervals for post preds (initial )
pp.pt <- data.frame(mn=numeric(n.preds),lci=numeric(n.preds),
                    uci=numeric(n.preds),dmu=numeric(n.preds),
                    year=numeric(n.preds),season=numeric(n.preds))
for(i in 1:n.preds){
    pp.pt[i,1] <- mean(pp[,i])
    pp.pt[i,2] <- quantile(pp[,i], prob=0.025)
    pp.pt[i,3] <- quantile(pp[,i], prob=0.975)
    pp.pt[i,4] <- ppidx[i,"dmu"]
    pp.pt[i,5] <- ppidx[i,"year"]
    pp.pt[i,6] <- ppidx[i,"season"]
}

y.idx <- data.frame(y=numeric(n.preds),dmu=numeric(n.preds),
                    year=numeric(n.preds),season=numeric(n.preds))
for(i in 1:n.preds){
  y.idx[i,"y"] <- y[as.numeric(ppidx[i,"season"]),
                  as.numeric(ppidx[i,"year"]),
                  as.numeric(ppidx[i,"dmu"])]
  y.idx[i,"dmu"] <- ppidx[i,"dmu"]
  y.idx[i,"year"] <- ppidx[i,"year"]
  y.idx[i,"season"] <- ppidx[i,"season"]
}

y.idx$lci <- rep(NA,nrow(y.idx))
y.idx$uci <- rep(NA, nrow(y.idx))
y.idx <- y.idx[,c("y",colnames(pp.pt)[-1])]
colnames(y.idx)[colnames(y.idx)=="y"] <- "val"
colnames(pp.pt)[colnames(pp.pt)=="mn"] <- "val"
y.idx$typ <- rep("data",nrow(y.idx))
pp.pt$typ <- rep("pred",nrow(pp.pt))
ypred <- rbind.data.frame(y.idx,pp.pt)

#points and CI's
# ggplot(df, aes(y=val,x=year,col=typ)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=lci,ymax=uci)) +
#   labs(y="Number of harvested deer",x="Year") +
#   # geom_point(data=df,mapping=aes(y=y,x=year,col="red")) +
#   facet_wrap(.~dmu) +
#   scale_color_manual(values=c("red","black"),
#                      label=c("Data","Prediction"),name="")+
#   theme(axis.text.x = element_text(angle=90,vjust=0.5))

# sites.sf <- st_read("./GIS Data/DMUs/dmus.shp",quiet=T)

ypred <- ypred %>% filter(year!=" 1")

g<-list()
  for(i in 1:length(unique(ypred$season))){
    df <- ypred %>% 
      filter(season==i)
    df$year <- factor(as.numeric(df$year),levels=1:dat$nPeriods)
    levels(df$year) <- dat$years
    df$typ <- factor(df$typ, levels=c("data","pred"))
    for(j in 1:length(unique(df$dmu))){
      df$dmu[df$dmu==j] <- sites.sf$Name[as.numeric(unique(df$dmu[df$dmu==j]))]
    }
    #lines and CI's
    df.p <- df %>% filter(typ=="pred")
    df.d <- df %>% filter(typ=="data")
    if(legend==T){
      g[[i]] <- ggplot(df.p, aes(y=val,x=year)) +
        geom_line(aes(group=1,color="Mean model estimate")) +
        geom_ribbon(aes(group=1,ymin=lci,ymax=uci,fill="95% credible interval"), alpha = 0.2) +
        geom_point(data=df.d,mapping=aes(y=val,x=year,col=typ)) +
        facet_wrap(.~dmu) +
        scale_colour_manual("",labels=c("Data","Mean model estimate"),
                            values=c("black","black"))+
        guides(colour = guide_legend(override.aes = list(pch = c(16, NA),lty=c(NA,1))))+
        scale_fill_manual("",values="blue")+
        theme(axis.text.x = element_text(angle=90,vjust=0.5)) +
        ggtitle(paste0(removal.names[i]," Harvest by DMU")) +
        labs(y="Deer Harvest",x="Year")
    }else {
      g[[i]]<- ggplot(df.p, aes(y=val,x=year)) +
        geom_line(aes(group=1,color="Mean model estimate")) +
        geom_ribbon(aes(group=1,ymin=lci,ymax=uci,fill="95% credible interval"), alpha = 0.2) +
        geom_point(data=df.d,mapping=aes(y=val,x=year,col=typ)) +
        facet_wrap(.~dmu) +
        scale_colour_manual("",labels=c("Data","Mean model estimate"),
                            values=c("black","black"))+
        guides(colour = "none",fill="none")+
        scale_fill_manual("",values="blue")+
        theme(axis.text.x = element_text(angle=90,vjust=0.5)) +
        ggtitle(paste0(removal.names[i]," Harvest by DMU")) +
        labs(y="Deer Harvest",x="Year")
    }
  }
  
  # g[[1]]
  # g[[2]]
  # 
  # gridExtra::grid.arrange(g[[1]],g[[2]])

if(nrow(y)==2){
  pp.plot <- gridExtra::grid.arrange(g[[1]],g[[2]],nrow=1)
} else if(nrow(y)==4){
  pp.plot <- gridExtra::grid.arrange(g[[1]],g[[2]],g[[3]],g[[4]],
                          nrow=2)
}
  list(pp.plot=pp.plot,df.p=df.p,data.pred.df=ypred)                
}

