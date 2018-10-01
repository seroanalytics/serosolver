setwd("~/net/home/coinflip")
source("scripts/cluster_setup.R")
source("scripts/run_test.R")
maxChain <- 10
runs <- read.csv("~/Documents/Fluscape/coinflip/testing_options.csv",stringsAsFactors=FALSE)
chainNos <- rep(1:maxChain,each=nrow(runs))
runs <- mefa:::rep.data.frame(runs,maxChain)
runs <- cbind(runs,"chainNo"=chainNos)
runs$runID <- "testing_gibbs"

iter <- 1000000
thin <- 1000
adapt_period <- 200000
adapt_freq <- 5000

indivs <- 100
n <- 50

pars <- c(4, 0.3, 1)
coin_probs <- runif(n,0,0.3)
samps <- seq(1,n, by=1)

coin_results <- sapply(coin_probs, function(x) sample(c(0,1),indivs,prob=c(1-x,x),replace=TRUE))
dat <- coin_toss_group(pars, coin_results)
dat <- measurement_error_group(pars,dat)

runs <- runs[,c("runID","chainNo","theta_proposal","lambda_proposal","Z_proposal","sampPropn","yearPropn")]
runs$theta_proposal <- as.character(runs$theta_proposal)
runs$lambda_proposal <- as.character(runs$lambda_proposal)
runs$Z_proposal <- as.character(runs$Z_proposal)
runs$sampPropn <- as.numeric(runs$sampPropn)
runs$yearPropn <- as.numeric(runs$yearPropn)
runs$chainNo <- as.numeric(runs$chainNo)
runs <- runs[runs$Z_proposal == "gibbs",]

jobs <- queuer::enqueue_bulk(obj, runs, "run_test",
                             swapPropn=0.5,adaptive=TRUE,pars=pars,
                             coin_probs=coin_probs,coin_results=coin_results,dat=dat,samps=samps,
                             iter=iter, thin=thin, adapt_period=adapt_period,adapt_freq=adapt_freq,plot=FALSE,
                             do_call=TRUE,timeout=0)
