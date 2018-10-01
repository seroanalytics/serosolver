######################
## JAMES HAY 13.08.2018 - jameshay218@gmail.com
## This script fits the serosolver antibody kinetics model to simulated HI titre data
## This particular script uses a gibbs proposal step to resample infection histories
## which integrates out the annual force of infection parameters.

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

filename <- "chains/sim_temp_test"

## How many individuals to simulate?
n_indiv <- 50
buckets <- 1 ## Set to 1 for annual model. Greater than 1 gives subannual (eg. buckets = 2 is infection period every half year)

## Read in parameter table to simulate from and change waning rate and alpha/beta if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)


samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)
yearMin <- 1990*buckets
yearMax <- 2015*buckets
ageMin <- 6*buckets
ageMax <- 75*buckets

## Read in and generate the antigenic map to read strain relationships from
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
fit_dat <- fit_dat[fit_dat$inf_years >= yearMin & fit_dat$inf_years <= yearMax,]
strainIsolationTimes <- unique(fit_dat$inf_years)

## Simulate data
dat <- simulate_data(parTab, 1, n_indiv, buckets, strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, ageMin,ageMax,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),
                     useSIR=TRUE, pInf = NULL, useSpline=FALSE)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets

titreDat <- dat[[1]]
titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

## If we want to use or save pre-simulated data
#write.table(titreDat,"~/net/home/serosolver/data/sim_1000_data.csv",sep=",",row.names=FALSE)
#write.table(infHist,"~/net/home/serosolver/data/sim_1000_infHist.csv",sep=",",row.names=FALSE)
#write.table(ages,"~/net/home/serosolver/data/sim_1000_ages.csv",sep=",",row.names=FALSE)
#write.table(AR,"~/net/home/serosolver/data/sim_1000_AR.csv",sep=",",row.names=FALSE)

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

## Generate starting locations for MCMC
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=50000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=10000,
              "save_block"=1000,"thin2"=100,"histSampleProb"=0.5,"switch_sample"=2, "burnin"=0, 
              "nInfs"=floor(ncol(infectionHistories)/2), "moveSize"=2*buckets, "histProposal"=6, "histOpt"=0,"n_par"=10, "swapPropn"=0.5)

## Run the MCMC using the inputs generated above
res <- run_MCMC(startTab, titreDat, mcmcPars, filename=filename,
                 create_post_func, NULL,NULL,version=1, 0.2, 
                 fit_dat, ages=ages, 
                 startInfHist=startInf,
                measurement_indices=NULL,
                temp=2)
beepr::beep(4)

#########################
## Processing outputs
#########################
## Density/trace plots
chain <- read.csv(res$chain_file)
chain <- chain[chain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]

## Plot inferred attack rates against true simulated attack rates
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]+10000),]
infChain <- infChain[infChain$sampno > 10000,]
n_alive <- sapply(1:length(strainIsolationTimes), function(x) length(ageMask[ageMask<=x]))
inf_prop <- colSums(infectionHistories)/n_alive
inf_prop <- data.frame(AR=inf_prop,year=strainIsolationTimes)
AR_p <- plot_attack_rates(infChain, titreDat, ages, seq(yearMin, yearMax, by=1)) + 
  geom_point(data=AR,aes(x=year,y=AR)) + 
  geom_point(data=inf_prop,aes(x=year,y=AR),col="purple") + scale_y_continuous(expand=c(0,0),limits=c(0,1))

## Density/trace plots on total number of infections
indivs <- sample(n_indiv, 10)
sampd <- sample(n_indiv,20)
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
n_strain <- max(infChain$j)
data.table::setkey(infChain, "j","sampno")
n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]

inf_chain_p <- ggplot(n_inf_chain) + geom_line(aes(x=sampno,y=V1)) + facet_wrap(~j)


n_indiv <- max(infChain$i)
data.table::setkey(infChain, "i","sampno")
i_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]

inf_chain_p_i <- ggplot(i_inf_chain) + geom_line(aes(x=sampno,y=V1)) + facet_wrap(~i)

## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
ps <- generate_cumulative_inf_plots(res$history_file, 0, 
                                    sampd, infHist, startInf,strainIsolationTimes, 100,ages,numberCol=4)


## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
plot_infection_histories(chain, infChain, titreDat, sample(1:n_indiv, 10), fit_dat, ages,parTab,100)

f <- create_post_func(parTab,titreDat, fit_dat, NULL, 99, ageMask, mu_indices=NULL,measurement_indices=NULL,temp=1000000)
tmpInf <- f(parTab$values, startInf,1,1,1,1,1,1)
samps <- 50000
tmpInf <- startInf
omg <- matrix(ncol=ncol(startInf)+1,nrow=nrow(startInf)*(samps))
omg <- matrix(ncol=nrow(startInf)+1,nrow=samps)
index <- 1
for(i in 1:samps){
  #print(i)
  tmp <- f(parTab$values, tmpInf, 1,1,1,1,3)
  tmpInf <- tmp
  tmp <- rowSums(tmp)
  tmp <- c(tmp,i)
  omg[i,] <- tmp
  #omg[index:(index+nrow(startInf)-1),] <- tmp
  #index <- index + nrow(startInf)
}
colnames(omg) <- c(1:nrow(startInf),"sampno")
omg <- as.data.frame(omg)
wow <- reshape2::melt(omg,id.vars="sampno")
ggplot(wow) + geom_histogram(aes(x=value),binwidth=2) + facet_wrap(~variable)

