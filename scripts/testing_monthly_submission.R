setwd("E:/James/Documents/Fluscape/serosolver")
devtools::load_all()

parTab <- read.csv("E:/James/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("E:/James/Documents/Fluscape/serosolver/data/fluscape_map.csv")
parTab[parTab$names == "wane","values"] <- 0.2

## How many individual to simulate/use? Leave as NULL if all individuals for real data
n_indiv <- 100


## Simulation options
samplingTimes <- (2007:2015)
nsamp <- 2

strainIsolationTimes <- unique(antigenicMap$inf_years)

## Ages between 5 and 80, censor 0% of titres
#dat <- simulate_data(parTab,1,n_indiv,strainIsolationTimes,samplingTimes, nsamp,antigenicMap, 0,0,5*12,80*12,simInfPars=c("mean"=0.01,"sd"=0.5,"bigMean"=0.08,"logSD"=1))
dat <- simulate_data(parTab,1,n_indiv,strainIsolationTimes,samplingTimes, nsamp,antigenicMap, 0,0,5,80)

## Extract simulation data
infectionHistories <- dat[["infectionHistories"]]
data <- dat[["data"]]
ages <- dat[["ages"]]
ages <- data.frame(individual=1:n_indiv,DOB=ages)
data$run <- 1
testedStrains <- seq(1970,2010,by=5)
#data <- data[data$virus %in% testedStrains,]

mcmcPars <- c("iterations"=50000,"popt"=0.44,"opt_freq"=1000,"thin"=1,"adaptive_period"=10000,
              "save_block"=100,"thin2"=10,"histSampleProb"=0.1,"switch_sample"=2, "burnin"=1000, 
              "nInfs"=2)

## For multivariate proposals
covMat <- diag(nrow(parTab))
scale <- 0.01
w <- 0.9
mvrPars <- list(covMat, scale, w)

## For univariate proposals
mvrPars <- NULL

f <- create_post_func(parTab,data,antigenicMap,NULL)

parTab[parTab$names == "error","fixed"] <- 1
startTab <- parTab
#f(parTab$values, infectionHistories)
for(i in 1:nrow(startTab)){
  if(startTab$fixed[i] == 0){
    startTab$values[i] <- runif(1,startTab$lower_bound[i],startTab$upper_bound[i])
  }
}


res <- lazymcmc::run_MCMC(startTab, data, mcmcPars, filename="test",create_post_func1, NULL, NULL, 0.2, 
                          antigenicMap=antigenicMap, infectionHistories=infectionHistories)

chain_lazy <- read.csv("test_univariate_chain.csv")
library(coda)
chain_lazy <- as.mcmc(chain_lazy[chain_lazy$sampno > 11000 & chain_lazy$sampno < 60002,])
#plot(coda::as.mcmc(chain))

res <- run_MCMC(startTab, data, mcmcPars, filename="test1",create_post_func, mvrPars, NULL, 0.2, antigenicMap, ages, startInfHist=infectionHistories)
chain_norm <- read.csv(res$chain_file)
chain_norm <- as.mcmc(chain_norm[chain_norm$sampno > 11000 & chain_norm$sampno < 60002,])

chains <- mcmc.list(chain_lazy,chain_norm)
plot(chains)
#plot(coda::as.mcmc(chain))
#plot(coda::as.mcmc(chain))

bestPars <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

startTab <- parTab
startTab$values <- bestPars
res <- run_MCMC(startTab, data, mcmcPars, filename="test1",create_post_func, mvrPars, NULL,0.2, antigenicMap, ages, startInfHist=infectionHistories)
chain <- read.csv(res$chain_file)
plot(coda::as.mcmc(chain))
