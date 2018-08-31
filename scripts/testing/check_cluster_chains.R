library(coda)
chains <- list.files(pattern="chain.csv")
omg <- NULL
for(chain in chains){
  tmp <- read.csv(chain)
  tmp <- tmp[tmp$sampno > 100000 & tmp$sampno < 500000, ]#c("mu","mu_short","tau","wane","sigma1","sigma2","error","lnlike")]
 
  #tmp <- tmp[tmp$sampno < 8000,c("mu","mu_short","tau","wane","sigma1","sigma2","error","lnlike")]
  omg[[chain]] <- as.mcmc(tmp)
}
omg <- as.mcmc.list(omg)
#omg <- as.mcmc.list(omg[c(1,2,5,6,8)])
#omg <- omg[c(-1,-6)]

pdf("tmp.pdf")
plot(omg)
dev.off()

library(ggplot2)
library(data.table)
library(plyr)
n_indiv <- 1000
indivs <- sample(n_indiv, 10)
sampd <- sample(n_indiv,20)
all_chain <- NULL
infHistFiles <- list.files(pattern="infectionHistories.csv")
for(i in 1:length(infHistFiles)){
  infChainFile <- infHistFiles[i]
  infChain <- data.table::fread(infChainFile)
   infChain <- infChain[infChain$sampno > 50000 & infChain$sampno < 100000,]
   n_strain <- max(infChain$j)
   data.table::setkey(infChain, "j","sampno")
   n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
   n_inf_chain <- cbind(n_inf_chain,"chain"=i)
   all_chain <- rbind(all_chain, n_inf_chain)
} 
#all_chain <- all_chain[all_chain$chain %in% c(1,2,5,6,8),]
all_chain$chain <- as.factor(all_chain$chain)
all_chain <- all_chain[!(all_chain$chain %in% c(1,6)),]
inf_chain_p <- ggplot(all_chain[all_chain$j %in% sample(unique(all_chain$j),50)]) + geom_line(aes(x=sampno,y=V1,col=chain)) + facet_wrap(~j)
devtools::load_all("~/Documents/Fluscape/serosolver")
titreFiles <- list.files(pattern="titreDat")
titreDat <- read.csv(titreFiles[1])
ageFiles <- list.files(pattern="ages")
ages <- read.csv(ageFiles[1])
mapFiles <- list.files(pattern="antigenic")
fit_dat <- read.csv(mapFiles[1])
strainIsolationTimes <- fit_dat$inf_years
AR_p <- plot_attack_rates_monthly(infChain, titreDat, ages, seq(1968*1,2015*1,by=1),ymax=1,buckets=1)
AR_p1 <- plot_attack_rates(infChain, titreDat, ages, seq(1968, 2015, by=1)) 
pdf("AR.pdf")
plot(AR_p)
dev.off()

pdf("AR1.pdf")
plot(AR_p1)
dev.off()

parTabFiles <- list.files(pattern="startTab")
parTab <- read.csv(parTabFiles[1])

chain <- read.csv(chains[1])
chain <- chain[chain$sampno > 500000 & chain$sampno < 1700000,]

clusters <- read.csv("~/Documents/Fluscape/serosolver/data/antigenic_maps/fonville_clusters.csv")
n_clusters <- length(unique(clusters$cluster1))
measurement_indices <- clusters$cluster1
infChain1 <- data.table::fread(infChainFile[1])
infChain1 <- infChain1[infChain1$sampno > 50000 & infChain1$sampno < 100000,]
pdf("fits.pdf")
plot_infection_histories(chain, infChain, titreDat, sample(1:100, 10), fit_dat, ages,parTab,nsamp=100, mu_indices=NULL,
                         measurement_indices=NULL)
dev.off()



infectionHistories <- infChain
dat <- titreDat
yearRange <- seq(1968*4,2015*4,by=1)
buckets <- 1
months <- 1:buckets
years <- range(floor(yearRange/buckets))
years <- years[1]:years[2]
labels <- c(sapply(years, function(x) paste0(months, "/",x)))
labels1 <- labels[1:length(yearRange)]
labels1 <- labels1[seq(1,length(labels1),by=buckets)]
yearBreak <- yearRange[seq(1,length(yearRange),by=buckets)]
n_alive <- sapply(yearRange, function(x) nrow(ages[ages$DOB <= x,]) )
##quantiles <- apply(tmp[,2:ncol(tmp)],2, function(x) quantile(x,c(0.025,0.5,0.975)))
data.table::setkey(infectionHistories, "sampno","j")
tmp <- infectionHistories[,list(V1=sum(x)),by=key(infectionHistories)]
##tmp <- ddply(infectionHistories, c("sampno","j"), function(x) sum(x$x))
quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025,0.5,0.975)))
colnames(quantiles) <- c("j","lower","median","upper")
quantiles[c("lower","median","upper")] <- quantiles[c("lower","median","upper")]/n_alive
##quantiles <- as.data.frame(t(quantiles))
quantiles$year <- yearRange[quantiles$j]
quantiles$taken <- quantiles$year %in% unique(dat$samples)

## Colour depending on whether or not titres were taken in each year
quantiles$taken <- ifelse(quantiles$taken,"Yes","No")
years <- rep(seq(1968,2015,by=1),each=4)
years <- years[1:nrow(quantiles)]
months <- rep(seq(1,4,by=1),length(seq(1968,2015,by=1)))
months <- months[1:nrow(quantiles)]
quantiles <- cbind(quantiles,"year"=years,"month"=months)
ggplot(quantiles) + geom_pointrange(aes(x=month,y=median,ymax=upper,ymin=lower,group=month))

