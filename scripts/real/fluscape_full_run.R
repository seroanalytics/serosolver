######################
## JAMES HAY 13.08.2018 - jameshay218@gmail.com
## This script fits the serosolver antibody kinetics model to the fluscape HI titre data
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
set.seed(1)

## How many individuals to fit to?
n_indiv <-50
buckets <- 1
## The general output filename
filename <- "chains/fluscape_quarter_test"

## Read in parameter table to simulate from and change waning rate and alpha/beta if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_base.csv",stringsAsFactors=FALSE)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
parTab <- parTab[parTab$names != "lambda",]


## Read in and generate the antigenic map to read strain relationships from
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, 1)
strainIsolationTimes <- unique(fit_dat$inf_years)


if("lambda" %in% parTab$names){
  ## Add rows for each lambda value to be inferred
  tmp <- parTab[parTab$names == "lambda",]
  for(i in 1:(length(strainIsolationTimes)-1)){
    parTab <- rbind(parTab, tmp)
  }
}
## Read in Fluscape data
fluscapeDat <- read.csv("data/real/combined_fluscape_data_12.csv",stringsAsFactors=FALSE)

## Remove individuals with NA for DOB
fluscapeDat <- fluscapeDat[complete.cases(fluscapeDat),]

## Take random subset of individuals
indivs <- sample(unique(fluscapeDat$individual),n_indiv)
indivs <- indivs[order(indivs)]
#indivs <- fluscapeAges[fluscapeAges$DOB >= 1990,"individual"]
#indivs <- unique(fluscapeAges$individual)
n_indiv <- length(indivs)

titreDat <- fluscapeDat[fluscapeDat$individual %in% indivs,]
titreDat$individual <- match(titreDat$individual, indivs)
titreDat <- titreDat[,c("individual", "samples", "virus", "titre", "run", "group","DOB")]

startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

## Multivariate proposals or univariate? Use univariate for now
covMat <- diag(nrow(parTab))
scale <- 0.00005
w <- 1
mvrPars <- list(covMat, scale, w)
mvrPars <- NULL
DOBs <- unique(titreDat[,c("individual","DOB")])[,2]
startInf <- setup_infection_histories_new(titreDat, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(DOBs, strainIsolationTimes)

#########################
## RUN MCMC
#########################
mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=50000,
              "save_block"=1000,"thin2"=100,"histSampleProb"=0.5,"switch_sample"=2, "burnin"=0, 
              "nInfs"=floor(ncol(startInf)/2), "moveSize"=2*buckets, "histProposal"=6, "histOpt"=0,"n_par"=10, "swapPropn"=0.5,
              "histSwitchProb"=0.2,"yearSwapPropn"=0.75)

res <- run_MCMC(startTab, titreDat=titreDat, antigenicMap=fit_dat,
                mcmcPars=c(save_block=1000, hist_switch_prob=0.5,year_swap_propn=0.5),
                mvrPars=NULL, filename=filename,
                create_posterior_func, PRIOR_FUNC=NULL,version=2, 0.2, 
                startInfHist=startInf,mu_indices=NULL,measurement_random_effects=FALSE,
                measurement_indices=NULL,
                temp=1)
beepr::beep(4)

#########################
## Processing outputs
#########################
chain <- read.csv(res$chain_file)
chain <- chain[chain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]

## Plot inferred attack rates
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
ageMask <- create_age_mask(DOBs, strainIsolationTimes)
strainMask <- create_strain_mask(titreDat,strainIsolationTimes)
masks <- data.frame(cbind(ageMask, strainMask))
n_alive <- sapply(seq(1,length(strainIsolationTimes)), function(x) nrow(masks[masks$ageMask <=x & masks$strainMask >= x,]))    

inf_prop <- colSums(startInf)/n_alive
inf_prop <- data.frame(AR=inf_prop,year=strainIsolationTimes)
AR_p <- plot_attack_rates(infChain, titreDat, ages, seq(min(strainIsolationTimes), max(strainIsolationTimes), by=1),n_alive=n_alive) + 
 scale_y_continuous(expand=c(0,0),limits=c(0,1))

## Density/trace plots on total number of infections
indivs <- sample(n_indiv, 10)
sampd <- sample(n_indiv,20)
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
n_strain <- max(infChain$j)
data.table::setkey(infChain, "j","sampno")
n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]

inf_chain_p <- ggplot(n_inf_chain[n_inf_chain$j < 10,]) + geom_line(aes(x=sampno,y=V1)) + facet_wrap(~j)

wow <- data.frame(n_inf_chain)
wow1 <- dcast(wow, sampno~j)
wow1[is.na(wow1)] <- 0

## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
plot_infection_histories(chain, infChain, titreDat, sample(1:n_indiv, 10), fit_dat, ages,parTab,100)

