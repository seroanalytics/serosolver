##### Sim recover vietnam like data
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
## How many individuals to simulate?
n_indiv <- 69

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- 1
## Possible sampling times
samplingTimes <- 2007:2012

## Antigenic map for cross reactivity parameters
fit_dat <- read.csv("~/Documents/Fluscape/serosolver/data/antigenicMap_AK.csv")

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968 & fit_dat$inf_years <= 2012,]
strainIsolationTimes <- unique(fit_dat$inf_years)

## Change alpha and beta to change proposal distribution
## Setting to c(1,1) gives uniform distribution on total number of infections
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(2,12)

## Simulate some fake data
hAR <- read.table("~/net/home/serosolver/data_LSA/HaNam_AR.txt",header=TRUE)
hAR <- hAR[,1]
dat <- simulate_data(parTab, 1, n_indiv, 1,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 6,75,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),
                     useSIR=FALSE, pInf = hAR, useSpline=FALSE)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets
#viruses <- seq(2000*buckets,2015*buckets,by=buckets)
#viruses <- seq(1968*buckets,2014*buckets,by=4)

titreDat <- dat[[1]]
#titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]


titreDat <- read.csv("data/vietnam_sim_dat_1.csv",stringsAsFactors=FALSE)
ages <- read.csv("data/vietnam_sim_ages_1.csv",stringsAsFactors=FALSE)
infectionHistories <- infHist <- read.csv("data/vietnam_sim_infHist_1.csv",stringsAsFactors=FALSE)
AR <- read.csv("data/vietnam_sim_AR_1.csv",stringsAsFactors=FALSE)
ages1 <- ages
ages1$DOB <- 1940

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

histProposal <- 4
version <- 1
## The general output filename
filename <- "chains/vietnam_original_sim"
## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,"adaptive_period"=50000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=1, "moveSize"=2, "histProposal"=histProposal, "histOpt"=0,"n_par"=10)

## Generate starting locations
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}
res <- run_MCMC(parTab, titreDat, mcmcPars, filename=filename,
                create_post_func, NULL,NULL,version=version, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)

infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= 50000,]
xs <- min(strainIsolationTimes):max(strainIsolationTimes)
n_alive <- sapply(strainIsolationTimes, function(x) nrow(ages[ages$DOB <= x,]))
arP <- plot_attack_rates(infChain, titreDat,ages,xs, n_alive)

colnames(AR) <- c("year","AR")
arP <- arP  + geom_point(data=AR,aes(x=year,y=AR),col="black",size=2,shape=1,stroke=1) + 
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) + 
  theme(text=element_text(family="Arial"),
        legend.position="none",
        axis.text.x=element_text(size=10,colour="black"),
        axis.text.y=element_text(size=10,colour="black"),
        axis.title.x=element_text(size=12,colour="black"),
        axis.title.y=element_text(size=12,colour="black")) +
  scale_color_manual(values=c("blue", "red", "black"))
arP



histProposal <- 1
version <- 4
## The general output filename
filename <- "chains/vietnam_correct_sim"
mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,"adaptive_period"=50000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=1, "moveSize"=2, "histProposal"=histProposal, "histOpt"=0,"n_par"=10)
## Add rows for each lambda value to be inferred
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_lambda.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- 1

tmp <- parTab[parTab$names == "lambda",]
for(i in 1:(length(strainIsolationTimes)-1)){
  parTab <- rbind(parTab, tmp)}

res <- run_MCMC(parTab, titreDat, mcmcPars, filename=filename,
                create_post_func, NULL,NULL,version=version, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)
