setwd("~/net/home/coinflip")
source("scripts/cluster_setup.R")
source("scripts/run_test.R")
maxChain <- 10
runs <- read.csv("~/Documents/Fluscape/coinflip/testing_options.csv",stringsAsFactors=FALSE)
chainNos <- rep(1:maxChain,each=nrow(runs))
runs <- mefa:::rep.data.frame(runs,maxChain)
runs <- cbind(runs,"chainNo"=chainNos)
runs$runID <- paste(runs$theta_proposal, runs$lambda_proposal, runs$Z_proposal,runs$sampPropn,runs$yearPropn,runs$chainNo,sep="_")

iter <- 50000
thin <- 100
adapt_period <- 10000
adapt_freq <- 1000

n_indiv <- 10
n <- 10

pars <- c(4, 0.3, 1)
coin_probs <- runif(n,0,0.3)

coin_results <- sapply(coin_probs, function(x) sample(c(0,1),indivs,prob=c(1-x,x),replace=TRUE))
dat <- coin_toss_group(pars, coin_results)
dat <- measurement_error_group(pars,dat)

runs <- runs[,c("runID","theta_proposal","lambda_proposal","Z_proposal","sampPropn","yearPropn")]

jobs <- queuer::enqueue_bulk(obj, runs, "run_test",
                             swapPropn=0.5,adaptive=TRUE,real_pars=pars,
                             coin_probs=coin_probs,coin_results=coin_results,dat=dat,samps=samps,
                             iter=iter, thin=thin, adapt_period=adapt_period,adapt_freq=adapt_freq,plot=FALSE,
                             do.call=TRUE,timeout=0)