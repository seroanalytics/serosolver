library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

n_indiv <- 10
version <- 1
histProposal <- 3
PRIOR <- infHistPrior

## Antigenic map fit
buckets <- 1
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets
antigenicMap$Strain <- virus_key[antigenicMap$Strain]
p1 <- ggplot(antigenicMap) + 
  geom_line(data=fit_dat,aes(x=x_coord,y=y_coord), col="red") +
  geom_point(data=antigenicMap,aes(x=X,y=Y)) + 
  geom_label(data=antigenicMap,aes(x=X+4,y=Y+0.25,label=Strain)) +
  theme_bw()

parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(2,18)
startTab <- parTab

fit_dat <- fit_dat[fit_dat$inf_years >= 1970*buckets,]
strainIsolationTimes <- unique(fit_dat$inf_years)
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 75*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0,"bigMean"=0.5,"logSD"=0),useSIR=TRUE)

titreDat <- dat[[1]]
titreDat <- titreDat[titreDat$virus %in% seq(1970*buckets,2010*buckets,by=buckets),]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

p <- plot_data(titreDat, infectionHistories, strainIsolationTimes, n_indiv, NULL)

startInf <- setup_infection_histories_new(titreDat, strainIsolationTimes, space=5,titre_cutoff=3)

optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
f1 <- create_post_func1(optimTab,titreDat,fit_dat,NULL,infectionHistories=startInf)
startPar <- parTab$values
startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=200))$optim$bestmem
startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])

mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=50000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=4, "moveSize"=5, "histProposal"=histProposal, "histOpt"=0)

## For univariate proposals

covMat <- diag(nrow(parTab))
scale <- 0.8
w <- 0.5
mvrPars <- list(covMat, scale, w)
f <- create_post_func(parTab,titreDat,fit_dat,version=version,PRIOR=PRIOR)
print(f(parTab$values, infHist))

mvrPars <- NULL

startTab$values <- startPar

#devtools::load_all()
res <- run_MCMC(startTab, titreDat, mcmcPars, filename="test1",
                create_post_func, mvrPars, PRIOR,version, 0.2, 
                fit_dat, ages=NULL, 
                startInfHist=startInf)


chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
plot(coda::as.mcmc(chain1))

infChain <- data.table::fread(res$history_file,data.table=FALSE)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
#infChain <- data.table::fread(infChainFile,data.table=FALSE)

n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
#beepr::beep(sound=4)
pdf(paste0(filename,"_infChain.pdf"))
plot(coda::as.mcmc(n_hist_chain))
dev.off()

n_hist_chain <- as.data.frame(n_hist_chain[,2:ncol(n_hist_chain)])
hist_chain <- reshape2::melt(n_hist_chain)


alpha <- parTab[parTab$names == "alpha","values"]
beta <- parTab[parTab$names == "beta","values"]

#y=dbb(0:(46*buckets), 46*buckets, alpha, beta)*dbinom(0:(46*buckets),46*buckets,0.5)
#y <- y/sum(y)
y=density_beta_binom(0:(46*buckets), 46*buckets, alpha, beta)
#y <- dbinom(0:(46*buckets),46,0.5)
dat <- data.frame(x=seq(0,46*buckets,by=1),y=y)

ggplot(hist_chain) + 
  geom_histogram(aes(x=value,y=..density..),binwidth=1,col="blue",fill="blue",alpha=0.1) + 
  facet_wrap(~variable) + 
  theme_bw() +
  geom_vline(data=totalInfsDat, aes(xintercept=y),col="red") +
  geom_line(data=dat,aes(x=x,y=y),col="red",size=0.5)

message(cat("Min ess theta: ", min(effectiveSize(chain1[,which(parTab$fixed==0)+1])),sep="\t"))
message(cat("Min ess infHist: ", min(effectiveSize(n_hist_chain)),sep="\t"))

#message(cat("Min gelman: ", min(effectiveSize(chain1[,which(parTab$fixed==0)+1])),sep="\t"))