library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
saveWD <- "~/Documents/Fluscape/serosolver_own/testing_redoc/"
filename <- "fluscape_sim_annual_measurement"

use_measurement_bias <- TRUE
## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1

## Antigenic map for cross reactivity parameters
#antigenic_map <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
#fit_dat <- generate_antigenic_map(antigenic_map, buckets)

fit_dat <- read.csv("data/antigenic_maps/created_maps/fonville_annual_continuous.csv")

## How many individuals to simulate?
n_indiv <- 100


## Read in parameter table to simulate from and change waning rate if necessary
par_tab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_base.csv",stringsAsFactors=FALSE)
par_tab[par_tab$names == "wane","values"] <- 0.8
par_tab[par_tab$names == "wane","values"] <- par_tab[par_tab$names == "wane","values"]/buckets

measurement_bias <- NULL
if(use_measurement_bias){
  measurement_bias <- rnorm(15,0,1)
  #measurement_bias[15] <- 0
  clusters <- read.csv("~/Documents/Fluscape/serosolver/data/antigenic_maps/fonville_clusters.csv")
  n_clusters <- length(unique(clusters$cluster1))
  measurement_indices <- clusters$cluster1
  measurement_indices <- rep(measurement_indices, each=buckets)

  for(i in 1:length(measurement_bias)){
    tmp <- data.frame(names="rho",values=1,fixed=0,steps=0.1,lower_bound=-10,upper_bound=10,lower_start=-2,upper_start=2, type=3)
    par_tab <- rbind(par_tab, tmp)
  }
  par_tab[par_tab$names == "rho","values"] <- measurement_bias
}

## Possible sampling times
sampling_times <- seq(2010*buckets, 2015*buckets, by=1)
#sampling_times <- 2007:2012
nsamps <- 2
repeats <- 1

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

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets & fit_dat$inf_years <= max(sampling_times),]
strain_isolation_times <- unique(fit_dat$inf_years)

#par_tab[par_tab$names %in% c("mu","mu_short","sigma1","sigma2"),"values"] <- c(2,2,0.3,0.1)
#par_tab[par_tab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strain_isolation_times),7,50)

## Simulate some fake data
## CHANGE PINF TO NULL IF WE WANT TO GENERATE NEW ATTACK RATES
sim_inf_pars=c("mean"=0.15/buckets,"sd"=0.5,"big_mean"=0.5,"log_sd"=1)
attack_rates <- simulate_attack_rates(strain_isolation_times, sim_inf_pars["mean"],sim_inf_pars["sd"],TRUE,sim_inf_pars["big_mean"])
dat <- simulate_data(par_tab, 1, n_indiv, buckets, strain_isolation_times,
                     sampling_times, nsamps, antigenic_map=fit_dat, titre_sensoring=0, 
                     age_min=10*buckets,age_max=75*buckets,
                     attack_rates=attack_rates,repeats=repeats, 
                     mu_indices=NULL,
                     measurement_indices = measurement_indices)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets

#viruses <- HaNam_viruses <- c(1968, 1972, 1976, 1982, 1989, 1992, 1993, 1994, 1995, 1996, 
#                   1997, 1999, 2000, 2002, 2003, 2004, 2005, 2007, 2008, 2009, 2010, 
#                   2011)*buckets

## If using HaNam data, need to filter the simulated data such that it matches the distribution of HaNam data
titre_dat <- dat[[1]]

## Create identifier for repeats
titre_dat$run <- NULL
titre_dat <- plyr::ddply(titre_dat,.(individual,virus,samples),function(x) cbind(x,"run"=1:nrow(x)))
titre_dat <- titre_dat[order(titre_dat$individual,titre_dat$run,titre_dat$samples,titre_dat$virus),]

## Merge with HaNam data to get same dimensions
#res <- read.csv("~/Documents/Fluscape/serosolver/data/real/vietnam_data.csv")
#wow <- dplyr::inner_join(res[,c("individual","samples","virus","run")], titre_dat[,c("individual","samples","virus","run")])
#titre_dat <- dplyr::left_join(wow, titre_dat)
#res <- read.csv("~/net/home/serosolver/data_LSA/HaNam_samples.csv")
#titre_dat <- merge(res[,c("individual","samples","virus","run")], titre_dat)
#titre_dat <- titre_dat[order(titre_dat$individual,titre_dat$run,titre_dat$samples,titre_dat$virus),]
titre_dat <- titre_dat[titre_dat$virus %in% viruses,]
infection_histories <- inf_hist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]
titre_dat <- merge(titre_dat,ages)

write.table(par_tab,paste0(saveWD,filename,"_pars_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(titre_dat,paste0(saveWD,filename,"_dat_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(inf_hist,paste0(saveWD,filename,"_inf_hist_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(ages,paste0(saveWD,filename,"_ages_",buckets,".csv"),row.names=FALSE,sep=",")
write.table(AR,paste0(saveWD,filename,"_AR_",buckets,".csv"),row.names=FALSE,sep=",")

