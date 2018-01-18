library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
monthly <- FALSE

filename <- "informative"
nInfs <- 10
filename <- paste0(filename,"_",nInfs)

if(!monthly){
  parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
  antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/fluscape_map.csv")
} else {
  parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_monthly.csv",stringsAsFactors=FALSE)
  antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/fluscape_map_quarter.csv")
}
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(2,10)



startTab <- parTab

## How many individual to simulate/use? Leave as NULL if all individuals for real data
n_indiv <- 10

## Simulation options
samplingTimes <- (2010:2015)
if(monthly) samplingTimes <- samplingTimes
nsamp <- 2

if(monthly){
  antigenicMap <- antigenicMap[antigenicMap$inf_years >= 1970*4,]
} else {
  antigenicMap <- antigenicMap[antigenicMap$inf_years >= 1970,]
}
strainIsolationTimes <- unique(antigenicMap$inf_years)

bucket <- 1
dat <- simulate_data(parTab, 1, n_indiv, bucket,strainIsolationTimes,
                     samplingTimes*bucket, 2, antigenicMap, 0, 0, 75*bucket,75*bucket,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=FALSE)
#plot(dat[[4]])
titreDat <- dat[[1]]
#titreDat <- titreDat[titreDat$virus %in% unique(antigenicMap$inf_years),]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

#parTab <- read.csv("~/net/home/serosolver/inputs/parTab.csv",stringsAsFactors = FALSE)
#antigenicMap <- read.csv("data/fluscape_map.csv")
#antigenicMap <- antigenicMap[antigenicMap$inf_years >= 1970,]
#titreDat <- read.csv("~/net/home/serosolver/data/sim_10_dat.csv",stringsAsFactors = FALSE)
#infHist <- read.csv("~/net/home/serosolver/data/sim_10_infHist.csv",stringsAsFactors = FALSE)
#ages <- read.csv("~/net/home/serosolver/data/sim_10_ages.csv",stringsAsFactors = FALSE)

#write.table(titreDat, paste0("~/net/home/serosolver/data/sim_",n_indiv,"_dat.csv"),row.names=FALSE, sep=",")
#write.table(infectionHistories, paste0("~/net/home/serosolver/data/sim_",n_indiv,"_infHist.csv"),row.names=FALSE, sep=",")
#write.table(ages, paste0("~/net/home/serosolver/data/sim_",n_indiv,"_ages.csv"),row.names=FALSE, sep=",")
#write.table(AR, paste0("~/net/home/serosolver/data/sim_",n_indiv,"_AR.csv"),row.names=FALSE, sep=",")

p <- plot_data(titreDat, infectionHistories, strainIsolationTimes, n_indiv, NULL)
#pdf(paste0("~/net/home/serosolver/data/sim_",n_indiv,"_data_p.pdf"))
#plot(p)
#dev.off()

startInf <- setup_infection_histories_new(titreDat, strainIsolationTimes, space=5,titre_cutoff=3)

optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
f1 <- create_post_func1(optimTab,titreDat,antigenicMap,NULL,infectionHistories=startInf)
startPar <- parTab$values
#startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=200))$optim$bestmem
#startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])

mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.234,"opt_freq"=2000,"thin"=10,"adaptive_period"=0,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=nInfs, "moveSize"=5)
              
            ## For univariate proposals

covMat <- diag(nrow(parTab))
scale <- 0.8
w <- 0.5
mvrPars <- list(covMat, scale, w)
f <- create_post_func(parTab,titreDat,antigenicMap,NULL)
mvrPars <- NULL

startTab$values <- startPar

res <- run_MCMC(startTab, titreDat, mcmcPars, filename="test1",
                create_post_func, mvrPars, NULL, 0.2, 
                antigenicMap, ages=NULL, 
                startInfHist=startInf)


chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
#plot(coda::as.mcmc(chain1))

infChain <- data.table::fread(res$history_file,data.table=FALSE)
#infChain <- data.table::fread(infChainFile,data.table=FALSE)
n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
#beepr::beep(sound=4)
pdf(paste0(filename,"_infChain.pdf"))
plot(coda::as.mcmc(n_hist_chain))
dev.off()


wow <- ddply(infChain, ~individual, function(x) colSums(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
wow <- reshape2::melt(wow, id.vars="individual")
wow$variable <- as.integer(wow$variable)

infHist1 <- as.data.frame(infHist)
infHist1$individual <- 1:n_indiv
colnames(infHist1) <- c(strainIsolationTimes,"individual")
infHist1 <- reshape2::melt(infHist1, id.vars="individual")
infHist1 <- infHist1[infHist1$value == 1,]

firstSamp <- infChain[infChain$sampno == 1,!(colnames(infChain) == "sampno")]
firstSamp <- reshape2::melt(firstSamp, id.vars="individual")
firstSamp <- firstSamp[firstSamp$value == 1,]

#histProfiles <- ddply(infChain, .(individual,sampno), function(x) cumsum(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
infChain <- infChain[infChain$sampno > mcmcPars["adaptive_period"],]
histProfiles <- t(apply(infChain, 1, function(x) cumsum(x[1:(ncol(infChain)-2)])))
histProfiles <- as.data.frame(cbind(histProfiles, "individual"=infChain[,c("individual")]))

histProfiles_lower <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.025))))
histProfiles_lower <- reshape2::melt(histProfiles_lower, id.vars="individual")
colnames(histProfiles_lower) <- c("individual","variable","lower")

histProfiles_upper <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.975))))
histProfiles_upper <- reshape2::melt(histProfiles_upper, id.vars="individual")
colnames(histProfiles_upper) <- c("individual","variable","upper")

histProfiles_median <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.5))))
histProfiles_median <- reshape2::melt(histProfiles_median, id.vars="individual")
colnames(histProfiles_median) <- c("individual","variable","median")

quantHist <- merge(histProfiles_lower, histProfiles_upper, by=c("individual","variable"))
quantHist <- merge(quantHist, histProfiles_median, by=c("individual","variable"))

realHistProfiles <- as.data.frame(t(apply(infHist, 1, cumsum)))
colnames(realHistProfiles) <- strainIsolationTimes
realHistProfiles$individual <- 1:n_indiv
realHistProfiles <- reshape2::melt(realHistProfiles,id.vars="individual")
  
startHistProfiles <- as.data.frame(t(apply(startInf, 1, cumsum)))
colnames(startHistProfiles) <- strainIsolationTimes
startHistProfiles$individual <- 1:n_indiv
startHistProfiles <- reshape2::melt(startHistProfiles,id.vars="individual")

n_indiv <- 10

p1 <- ggplot(quantHist[quantHist$individual %in% seq(1,n_indiv,by=1),]) + 
  geom_line(aes(x=as.integer(variable),y=median)) + 
  geom_ribbon(aes(x=as.integer(variable),ymin=lower,ymax=upper), alpha=0.2) + 
  geom_line(data=realHistProfiles[realHistProfiles$individual %in% seq(1,n_indiv,by=1),],aes(x=as.integer(variable),y=value), col="blue") +
  geom_line(data=startHistProfiles[startHistProfiles$individual %in% seq(1,n_indiv,by=1),],aes(x=as.integer(variable),y=value), col="red") +
  #scale_y_continuous(limits=c(0,10)) +
  facet_wrap(~individual, scales="free_y")

p2 <- ggplot(wow) + geom_line(aes(x=variable,y=value)) + 
    geom_vline(data=infHist1,aes(xintercept=as.integer(variable)),col="red") + 
  #geom_vline(data=firstSamp, aes(xintercept=as.integer(variable)),col="blue") +
    facet_wrap(~individual)

pdf(paste0(filename,"_infHistP.pdf"))
plot(p1)
dev.off()

pdf(paste0(filename,"_infHistP2.pdf"))
plot(p2)
dev.off()

#pdf("original_proposal_prior.pdf")
#plot(coda::as.mcmc(n_hist_chain))
#dev.off()

#pdf("original_cum_inf.pdf")
#plot(p1)
#dev.off()

#pdf("original_each_inf.pdf")
#plot(p2)
#dev.off()
