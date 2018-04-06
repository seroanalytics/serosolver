library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
saveWD <- "~/net/home/serosolver/data_LSA/"
filename <- "vietnam_sim"

## How many individuals to simulate?
n_indiv <- 69

## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- 1
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets


## Possible sampling times
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)
samplingTimes <- 2007:2012
nsamps <- 6
repeats <- 6

############################
## SIMULATE ATTACK RATES
############################
## If using HaNam sim, use ARs from Adam's code
hAR <- read.table("~/net/home/serosolver/data_LSA/HaNam_AR.txt",header=TRUE)
hAR <- hAR[,1]
#hAR <- c(0.654, 0.965, 0.392, 0.215, 0.131, 0.0734, 0.107, 0.139, 0.13, 
#         0.0557, 0.0173, 0.0086, 0.0196, 0.00401, 0.00403, 0.0499, 0.222, 
#         0.799, 0.0193, 0.00724, 0.0112, 0.00443, 0.00369, 0.0282, 0.491, 
#         0.0987, 0.0233, 0.00247, 0.00551, 0.00812, 0.284, 0.155, 0.0311, 
#         0.637, 0.932, 0.349, 0.0241, 0.00303, 0.00376, 0.0303, 0.157, 
#         0.228, 0.104, 0.328, 0.306, 0.431, 0.4, 0.377)
#hAR <- NULL
#sdPar=0.5
#hAR <- rlnorm(48, meanlog=log(0.15)-(sdPar/2)^2/2,sdlog=sdPar/2)
#hAR[c(1,3)] <- rlnorm(2, meanlog=log(0.5)-(sdPar/2)^2/2,sdlog=sdPar/2)
#hAR[c(2,4)] <- rlnorm(2, meanlog=log(0.05)-(sdPar/2)^2/2,sdlog=sdPar/2)


## Antigenic map for cross reactivity parameters
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
fit_dat <- read.csv("data/antigenicMap_AK.csv")

## Rename circulation years based on isolation time
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets & fit_dat$inf_years <= max(samplingTimes),]
strainIsolationTimes <- unique(fit_dat$inf_years)

#parTab[parTab$names %in% c("mu","mu_short","sigma1","sigma2"),"values"] <- c(2,2,0.3,0.1)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)

## Simulate some fake data
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, nsamps, antigenicMap=fit_dat, 0, 0, 10*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=TRUE,pInf=hAR,
                     useSpline=FALSE,
                     repeats=repeats)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets

HaNam_viruses <- c(1968, 1972, 1976, 1982, 1989, 1992, 1993, 1994, 1995, 1996, 
                   1997, 1999, 2000, 2002, 2003, 2004, 2005, 2007, 2008, 2009, 2010, 
                   2011)

## If using HaNam data, need to filter the simulated data such that it matches the distribution of HaNam data
titreDat <- dat[[1]]

## Create identifier for repeats
titreDat$run <- NULL

titreDat <- plyr::ddply(titreDat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x)))
#for(indiv in unique(titreDat$individual)){
#  print(indiv)
#  tmp1 <- titreDat[titreDat$individual == indiv,]
#  for(sample in unique(tmp1$samples)){
#    tmp2 <- tmp1[tmp1$samples == sample,]
#    for(strain in unique(tmp2$virus)){
#      tmp3 <- tmp2[tmp2$virus == strain,"virus"]
#      tmp_seq <- seq(1,length(tmp3),by=1)
#      titreDat[titreDat$individual == indiv & titreDat$samples == sample & titreDat$virus == strain,"run"] <- tmp_seq
#    }
#  }
#}

titreDat <- titreDat[order(titreDat$individual,titreDat$run,titreDat$samples,titreDat$virus),]

## Merge with HaNam data to get same dimensions
res <- read.csv("~/Documents/Fluscape/serosolver/data/vietnam_data.csv")
wow <- dplyr::inner_join(res[,c("individual","samples","virus","run")], titreDat[,c("individual","samples","virus","run")])
titreDat <- dplyr::left_join(wow, titreDat)
#res <- read.csv("~/net/home/serosolver/data_LSA/HaNam_samples.csv")
#titreDat <- merge(res[,c("individual","samples","virus","run")], titreDat)
#titreDat <- titreDat[order(titreDat$individual,titreDat$run,titreDat$samples,titreDat$virus),]
#titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

saveWD <- "~/Documents/Fluscape/serosolver/data/"
write.table(parTab,paste0(saveWD,filename,"_pars_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(titreDat,paste0(saveWD,filename,"_dat_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(infHist,paste0(saveWD,filename,"_infHist_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(ages,paste0(saveWD,filename,"_ages_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(AR,paste0(saveWD,filename,"_AR_",buckets,".csv"),row.names=FALSE,sep=",")
