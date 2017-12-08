library(doMC)
library(foreach)
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
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/fluscape_map.csv")

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
    dat <- simulate_data(parTab,1,n_indiv,strainIsolationTimes,
                         samplingTimes, nsamp,antigenicMap, 0,0,5,80)

    ## Extract simulation data
    infectionHistories <- dat[["infectionHistories"]]
    data <- dat[["data"]]
    ages <- dat[["ages"]]
    ages <- data.frame(individual=1:n_indiv,DOB=ages)
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
mcmcPars <- c("iterations"=500000,"popt"=0.44,"opt_freq"=1000,"thin"=10,"adaptive_period"=100000,
              "save_block"=100,"thin2"=1000,"histSampleProb"=0.1,"switch_sample"=2, "burnin"=100000)

## For multivariate proposals
covMat <- diag(nrow(parTab))
scale <- 0.01
w <- 0.9
mvrPars <- list(covMat, scale, w)

## For univariate proposals
mvrPars <- NULL

## Generate random starting points and run
res <- foreach(x =filenames) %dopar% {
    startTab <- parTab
    for(i in 1:nrow(startTab)){
        if(startTab$fixed[i] == 0){
            startTab$values[i] <- runif(1,startTab$lower_bound[i],startTab$upper_bound[i])
        }
    }
    run_MCMC(startTab, data, mcmcPars, filename="test",create_post_func, NULL, mvrPars, 0.2, antigenicMap, ages, startInfHist=NULL)
}

res <- run_MCMC(startTab, data, mcmcPars, filename="test1",create_post_func, NULL, mvrPars, 0.2, antigenicMap, ages, startInfHist=NULL)

output <- res[[1]]
generate_all_plots(outputDir, mcmcPars["adaptive_period"], output$chain_file, output$history_file,
                   data, antigenicMap, parTab, ages, 10, 1000, "testing")


generate_all_plots(outputDir, mcmcPars["adaptive_period"], "test_chain.csv", "test_infectionHistories.csv",
                   data, antigenicMap, parTab, ages, 10, 100, "testing")
chain <- read.csv("test_chain.csv")
infectionHistories <- read.csv("test_infectionHistories.csv")
