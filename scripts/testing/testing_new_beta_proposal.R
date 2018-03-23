library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

## How many individuals to simulate?
n_indiv <-100

## Which infection history proposal version to use?
describe_proposals()
histProposal <- 1

## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1

## The general output filename
filename <- "chains/test_lambda1"

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_lambda.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- 1
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets
parTab[parTab$names == "sigma1","values"] <- parTab[parTab$names == "sigma1","values"]#*buckets
parTab[parTab$names == "sigma2","values"] <- parTab[parTab$names == "sigma2","values"]#*buckets
#parTab[parTab$names == "error","fixed"] <- 1

## Possible sampling times
samplingTimes <- seq(1968*buckets, 2015*buckets, by=1)

## Antigenic map for cross reactivity parameters
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)

## Rename circulation years based on isolation time
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets
antigenicMap$Strain <- virus_key[antigenicMap$Strain]

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets & fit_dat$inf_years <= 2015*buckets,]
strainIsolationTimes <- unique(fit_dat$inf_years)

## Add rows for each lambda value to be inferred
tmp <- parTab[parTab$names == "lambda",]
for(i in 1:(length(strainIsolationTimes)-1)){
  parTab <- rbind(parTab, tmp)
}

## Change alpha and beta to change proposal distribution
## Setting to c(1,1) gives uniform distribution on total number of infections
#parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)

## Simulate some fake data
#strainIsolationTimes <- 1968:2015
lambdas <- runif(length(strainIsolationTimes),0.1/buckets,0.5/buckets)
lambdas[1] <- 0.6/buckets
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 6*buckets,80*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),
                     useSIR=TRUE, pInf = NULL, useSpline=FALSE)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets
#viruses <- seq(1968*buckets,2015*buckets,by=buckets)
titreDat <- dat[[1]]
titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

#write.table(titreDat,paste0("~/net/home/serosolver/data/sim_lambda_1000_data_",buckets,"_.csv"),sep=",",row.names=FALSE)
#write.table(infHist,paste0("~/net/home/serosolver/data/sim_lambda_1000_infHist_",buckets,"_.csv"),sep=",",row.names=FALSE)
#write.table(ages,paste0("~/net/home/serosolver/data/sim_lambda_1000_ages_",buckets,"_.csv"),sep=",",row.names=FALSE)
#write.table(AR,paste0("~/net/home/serosolver/data/sim_lambda_1000_AR_",buckets,"_.csv"),sep=",",row.names=FALSE)
#write.table(parTab,paste0("~/net/home/serosolver/data/parTab_callibration_",buckets,"_.csv"),sep=",",row.names=FALSE)

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

## Generate starting locations
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=50000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,"adaptive_period"=20000,
              "save_block"=1000,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=20000, 
              "nInfs"=3, "moveSize"=5, "histOpt"=0,"nYears"=round(length(lambdas)*0.05))
covMat <- diag(nrow(parTab[parTab$block != 2,]))
scale <- 0.5
w <- 1
mvrPars <- list(covMat, scale, w)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)
f <- create_post_func(parTab=parTab,data=titreDat,antigenicMap=fit_dat,PRIOR_FUNC = NULL,version=4, ageMask=ageMask)
f_dat <- create_post_func(parTab=parTab,data=titreDat,antigenicMap=fit_dat,PRIOR_FUNC = NULL,version=1, ageMask=ageMask)
f_lambda <- create_post_func(parTab=parTab,data=titreDat,antigenicMap=fit_dat,PRIOR_FUNC = NULL,version=3, ageMask=ageMask)
parTab[parTab$names == "lambda","values"] <- lambdas
sum(f(parTab$values, infHist))
sum(f_dat(parTab$values, infHist))
sum(f_lambda(parTab$values, infHist))
parTab[parTab$names %in% c("alpha","beta"),"values"] <- startTab[startTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
lambdas_start <- NULL
for(year in 1:ncol(startInf)){ lambdas_start[year] <- (sum(startInf[,year])/sum(ages$DOB <= strainIsolationTimes[year]))}
#startTab[startTab$names == "lambda","values"] <- lambdas_start
parTab[parTab$names == "lambda","block"] <- 2
#parTab[parTab$names == "lambda","values"] <- 0.8
#startInf <- matrix(sample(c(0,1),n_indiv*length(strainIsolationTimes),replace=TRUE),nrow=n_indiv)
#for(i in 1:length(ageMask)){
#  if(ageMask[i] > 1){
#    startInf[i,1:(ageMask[i])] <- 0
#  }
#}
parTab[parTab$names == "lambda","values"] <- lambdas
#startTab[1:10,] <- parTab[1:10,]
parTab[parTab$names %in% c("alpha","beta"),"values"] <- startTab[startTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
Rprof(tmp <- tempfile())
res <- run_MCMC_lambda(startTab, titreDat, mcmcPars, filename=filename,
                OPT_TUNING=0.2, 
                antigenicMap=fit_dat, ages=ages, 
                startInfHist=startInf,
                ver=4,mvrPars=NULL,
                block_weights=c(2,2,1))
Rprof()
summaryRprof(tmp)
library(proftools)
plotProfileCallGraph(readProfileData(tmp),score = "total")

chain1 <-read.csv(res$chain_file)
chain1 <- chain1[seq(1,nrow(chain1),by=100),]
#chain1 <- chain1[chain1$sampno > 40000,]
liks_right <- NULL
for(i in 1:ncol(infHist)){
  liks_right[i] <- dbinom(sum(infHist[,i]),size=nrow(infHist),prob=lambdas[i],log=FALSE)
}
tmp <- chain1[,12:(length(lambdas)+11)]

liks <- matrix(nrow=nrow(tmp),ncol=ncol(infHist))
for(j in 1:nrow(tmp)){
  liksA <- NULL
  for(i in 1:ncol(infHist)){
    liksA[i] <- dbinom(sum(infHist[,i]),size=nrow(infHist),prob=tmp[j,i],log=FALSE)
  }
  liks[j,] <- liksA
}
wow <- melt(liks)
colnames(wow) <- c("sampno","variable","prob")
wow$variable <- colnames(tmp)[wow$variable]
tmp1 <- chain1[,c(1,12:(length(lambdas)+11))]
wow1 <- melt(tmp1,id.vars="sampno")
wow$sampno <- unique(tmp1$sampno)[wow$sampno]
omg <- merge(wow, wow1, by=c("sampno","variable"))
likDat <- data.frame(variable=1:length(liks_right),y=lambdas,y1=liks_right)
likDat$variable <- colnames(tmp)[likDat$variable]
p1 <- ggplot(omg) + geom_line(aes(x=value,y=prob)) + geom_vline(data=likDat,aes(xintercept=y),col="red") + facet_wrap(~variable)
p2 <- ggplot(wow) + geom_histogram(aes(x=prob)) + geom_vline(data=likDat,aes(xintercept=y1),col="red") + facet_wrap(~variable)

plot(coda::as.mcmc(chain1[,c("prob_full","prob_dat","prob_lambda")]))
chain1 <- read.csv(res1$chain_file)

ess <- effectiveSize(chain)
ess1 <- effectiveSize(chain1)


mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,"adaptive_period"=50000,
              "save_block"=1000,"thin2"=100,"histSampleProb"=0.1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=1, "moveSize"=2, "histProposal"=1, "histOpt"=0,"n_par"=1)
res1 <- run_MCMC(startTab, titreDat, mcmcPars, filename="chains/original_prop",
                 create_post_func, NULL,NULL,version=4, 0.2, 
                 fit_dat, ages=ages, 
                 startInfHist=startInf)

pars <- c("mu","mu_short","wane","sampno","tau","sigma1","sigma2","error","lambda","lambda.1")
melted_chain <- melt(chain1[chain1$sampno > 20000,pars],id.vars="sampno")
labels <- parTab[parTab$names %in% pars,c("names","values","lower_bound","upper_bound")]
colnames(labels) <- c("variable","values","lower_bound","upper_bound")
omg <- melt(labels[,c("variable","lower_bound","upper_bound")], id.vars="variable")
omg <- omg[,c(1,3)]
ggplot(melted_chain) + geom_density(aes(x=value,y=..scaled..)) + 
  facet_wrap(~variable,scales="free_x") + theme_bw()+
  geom_blank(data=omg,aes(x=value))+
  geom_vline(data=labels,aes(xintercept=values),col="red")
  