######################
## KYLIE AINSLIE 21.11.2018 - ainslie.kylie@gmail.com
## Minimal script for running serosolver using simulated data

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

# Set working directory and load code from serosolver package
# setwd("Volumes/Home/kainslie/serosolver") # Mac path
  setwd("Q:/serosolver") # PC path
  devtools::load_all()

###################################      
# Simulate H3N2 data #
###################################
# How many individuals to simulate?
  n_indiv <- 100

  buckets <- 1 ## Set to 1 for annual model. Greater than 1 gives subannual (eg. buckets = 2 is infection period every half year)

# Read in parameter table to simulate from and change waning rate and alpha/beta if necessary
  parTab <- read.csv("Q:/serosolver_testing/parTab_base.csv",stringsAsFactors=FALSE)
  parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
  parTab[parTab$names == "wane","values"] <- 1
  parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets
  parTab <- parTab[parTab$names != "phi",]

  samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)
  yearMin <- 1968*buckets
  yearMax <- 2015*buckets
  ageMin <- 6*buckets
  ageMax <- 75*buckets

# Read in and generate the antigenic map to read strain relationships from
  antigenicMap <- read.csv("Q:/serosolver_testing/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
  fit_dat <- generate_antigenic_map(antigenicMap, buckets)
  fit_dat <- fit_dat[fit_dat$inf_years >= yearMin & fit_dat$inf_years <= yearMax,]
  strainIsolationTimes <- unique(fit_dat$inf_years)

  simInfPars=c("mean"=0.15/buckets,"sd"=0.5,"bigMean"=0.5/buckets/2,"logSD"=1)
# use SIR
# attackRates  <- simulate_ars_buckets(strainIsolationTimes, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
# attackRates  <- simulate_attack_rates(strainIsolationTimes, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])
# attackRates  <- simulate_ars_spline(strainIsolationTimes, buckets, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"], knots,theta)
  attackRates <- simulate_attack_rates(strainIsolationTimes, simInfPars["mean"],simInfPars["sd"],TRUE,simInfPars["bigMean"])

# Simulate data
  dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, ageMin*buckets,ageMax*buckets,
                     attackRates)


###################################      
# Running serosolver on H3N2 data #
###################################

# The general output filename
  filename <- "test"

# add DOB column
  titreDat <- merge(dat$data,dat$ages,by='individual')
  
# Uncomment if you want to take a random subset of individuals
# indivs <- sample(unique(titreDat$individual),n_indiv)
# indivs <- indivs[order(indivs)]
  indivs <- unique(titreDat$individual) #all individuals

# Subset data for indivs (if not using all individuals)
#  titreDat <- titreDat[titreDat$individual %in% indivs,]
#  titreDat$individual <- match(titreDat$individual, indivs)

startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

# Multivariate proposals or univariate? Use univariate for now
  covMat <- diag(nrow(parTab))
  scale <- 0.00005
  w <- 1
  mvrPars <- list(covMat, scale, w)
  mvrPars <- NULL

# unique rows for each individual
  unique_indiv <- titreDat[!duplicated(titreDat$individual),]
  startInf <- setup_infection_histories_new(titreDat, strainIsolationTimes, space=10,titre_cutoff=2)
  ageMask <- create_age_mask(unique_indiv$DOB, strainIsolationTimes)

############
# RUN MCMC #
############

# default mcmcPars (note: don't need to pass to run_MCMC unless you want non-default options)
mcmcPars <- c("iterations"=50000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,"adaptive_period"=10000, 
              "save_block"=100, "thin2"=10,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, "inf_propn"=1, 
              "moveSize"=5,"histOpt"=1,"swapPropn"=0.5,"hist_switch_prob"=0,"year_swap_propn"=1)

res <- run_MCMC(parTab = startTab, titreDat = titreDat,antigenicMap = fit_dat,startInfHist=startInf, 
                filename = filename, CREATE_POSTERIOR_FUNC = create_posterior_func, version = 2)

beepr::beep(4)

######################
# Processing outputs #
######################

# Density/trace plots
  chain1 <- read.csv(res$chain_file)
  chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
  plot(coda::as.mcmc(chain1))
  
# Plot inferred attack rates
  infChain <- data.table::fread(res$history_file,data.table=FALSE)
  infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
  infChain <- setDT(infChain) # change to data.table
  xs <- min(strainIsolationTimes):max(strainIsolationTimes)
  AR_plot <- plot_attack_rates(infectionHistories = infChain, dat = titreDat, ages = titreDat[,c('individual', 'DOB')],
                               yearRange = xs, buckets=12)
# Note: for monthly data (buckets=12) use plot_attack_rates_monthly()

# Plot infection histories
  IH_plot <- plot_infection_histories(chain = chain1, infectionHistories = infChain, titreDat = titreDat, individuals = c(1:5),
                                      antigenicMap = fit_dat, parTab = startTab1)
