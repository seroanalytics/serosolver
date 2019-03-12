library(doMC)
library(foreach)
library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
#library(serosolver)

## Where to save final plots
outputDir <- "outputs"

## Simulate or real?
SIM <- TRUE

## How many individual to simulate/use? Leave as NULL if all individuals for real data
n_indiv <- 10

## CHANGE FOR LOCAL FILE SYSTEM
## Important input parameters and antigenic map
par_tab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_base.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/antigenicMap_AK.csv")

## Simulation options
samplingTimes <- 2007:2015
nsamp <- 2

## CHANGE FOR LOCAL FILE SYSTEM
## Make up filenames to save output to
filenames <- c("chains/output","chains/output","chains/output")

## We'll be parallelising a few chains
registerDoMC(cores=4)

if(SIM){
    ## Extract possible infection times
    strainIsolationTimes <- unique(antigenicMap$inf_years)

    ## Ages between 5 and 80, censor 0% of titres
    dat <- simulate_data(par_tab, 1, n_indiv, buckets,strainIsolationTimes,
                         samplingTimes, 2, antigenicMap=antigenicMap, 0, 0, 10,75,
                         simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=FALSE)

    ## Extract simulation data
    infectionHistories <- dat[["infectionHistories"]]
    data <- dat[["data"]]
    ages <- dat[["ages"]]
} else {
    ## CHANGE FOR LOCAL FILE SYSTEM
    data <- read.csv("data/fluscape_data.csv",stringsAsFactors=FALSE)
    ages <- read.csv("data/fluscape_ages.csv")

    if(!is.null(n_indiv)){
        indivs <- sample(unique(data$individual),n_indiv)
        data <- data[data$individual %in% indivs,]
        ages <- ages[ages$individual %in% indivs,]
    }
}


## MCMC parameter inputs
mcmcPars <- c("iterations"=20000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=1000,"thin"=1,"adaptive_period"=10000,
"save_block"=100,"thin2"=1,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
"nInfs"=4, "moveSize"=5, "histProposal"=3, "histOpt"=1)

## For multivariate proposals
covMat <- diag(nrow(par_tab))
scale <- 0.01
w <- 0.9
mvrPars <- list(covMat, scale, w)

## For univariate proposals
mvrPars <- NULL

## Generate random starting points and run
res <- foreach(x =filenames) %dopar% {
    startTab <- par_tab
    for(i in 1:nrow(startTab)){
        if(startTab$fixed[i] == 0){
            startTab$values[i] <- runif(1,startTab$lower_bound[i],startTab$upper_bound[i])
        }
    }
    run_MCMC(startTab, data, mcmcPars, filename="test",create_post_func, NULL, version = 1, mvrPars, 0.2, 
             antigenicMap, ages, startInfHist=NULL)
}
startTab <- par_tab
res <- run_MCMC(startTab, data, mcmcPars, filename="test1",create_post_func, NULL, version = 1, mvrPars, 0.2, 
                antigenicMap, ages, startInfHist=NULL)

output <- res
generate_all_plots(outputDir, mcmcPars["adaptive_period"], output$chain_file, output$history_file,
                   data, antigenicMap, par_tab, ages, 10, 100, "testing")
