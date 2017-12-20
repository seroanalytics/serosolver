library(ggplot2)
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()


parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_monthly.csv",stringsAsFactors=FALSE)
antigenicMap <- read.csv("~/Documents/Fluscape/serosolver/data/fluscape_map_monthly.csv")

## How many individual to simulate/use? Leave as NULL if all individuals for real data
n_indiv <- 1

## Simulation options
samplingTimes <- (2010:2015)*12
nsamp <- 5

antigenicMap <- antigenicMap[antigenicMap$inf_years >= 1970*12,]
strainIsolationTimes <- unique(antigenicMap$inf_years)

## Ages between 5 and 80, censor 0% of titres
## Create antigenic map for short and long term boosting
pars <- parTab$values
names(pars) <- parTab$names

antigenicMap1 <- outputdmatrix.fromcoord(antigenicMap[,c("x_coord","y_coord")])
antigenicMapLong <- 1 - pars["sigma1"]*c(antigenicMap1)
antigenicMapShort <- 1 - pars["sigma2"]*c(antigenicMap1)
antigenicMapLong[antigenicMapLong < 0] <- 0
antigenicMapShort[antigenicMapShort < 0] <- 0


infHist <- matrix(0, nrow=n_indiv, ncol=length(strainIsolationTimes))
for(i in 1:nrow(infHist)){
  infections <- sample(1:ncol(infHist),rpois(1, 5))
  infHist[i,infections] <- 1
}

data <- simulate_group(n_indiv, pars, infHist, strainIsolationTimes, samplingTimes, nsamp, antigenicMapLong, antigenicMapShort)
data$run <- 1
#testedStrains <- seq(1970*12,2010*12,by=24)
data <- data[data$virus %in% testedStrains,]
dat_plot1 <- ggplot(data[data$individual <= 5,]) + geom_point(aes(x=as.integer(virus),y=titre)) + facet_grid(individual~samples)

startInf <- matrix(0, nrow=n_indiv, ncol=length(strainIsolationTimes))
for(i in 1:nrow(infHist)){
  infections <- sample(1:ncol(infHist),rpois(1, 3))
  startInf[i,infections] <- 1
}


mcmcPars <- c("iterations"=100000,"popt"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=100000,
              "save_block"=100,"thin2"=100,"histSampleProb"=0.05,"switch_sample"=2, "burnin"=0, 
              "nInfs"=3)
              
            ## For univariate proposals

covMat <- diag(nrow(parTab))
scale <- 0.01
w <- 0.99
mvrPars <- list(covMat, scale, w)
f <- create_post_func(parTab,data,antigenicMap,NULL)
#mvrPars <- NULL

set.seed(2)
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab$fixed[i] == 0){
    startTab$values[i] <- runif(1,startTab$lower_bound[i],startTab$upper_bound[i])
  }
}
#startTab[startTab$names == "mu","values"] <- 2
#startTab[startTab$names == "mu_short","values"] <- 2
ages <- data.frame(individual=1:n_indiv, DOB=1970*12)

#res <- lazymcmc::run_MCMC(startTab, data, mcmcPars, filename="test",create_post_func, NULL, NULL, 0.2, 
#                          antigenicMap=antigenicMap, infectionHistories=infectionHistories)
res <- run_MCMC(startTab, data, mcmcPars, filename="test1",
                create_post_func, mvrPars, NULL, 0.2, 
                antigenicMap, ages, 
                startInfHist=startInf)

#chain <- read.csv(res$chain_file)
#bestI <- floor((which.max(chain$lnlike) / 10))*10 + 1
#bestPars <- get_best_pars(chain)
#chain <- chain[chain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),2:(ncol(chain)-1)]
#covMat <- cov(chain)
#mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

#startTab <- parTab
#startTab$values <- bestPars

#infChain <- data.table::fread(res$history_file,data.table=FALSE)
#bestInf <- as.matrix(infChain[infChain$sampno == bestI, 1:(ncol(infChain)-2)])
#mcmcPars <- c("iterations"=100000,"popt"=0.234,"opt_freq"=2000,"thin"=1,"adaptive_period"=20000,
#              "save_block"=100,"thin2"=10,"histSampleProb"=1,"switch_sample"=2, "burnin"=1000, 
#              "nInfs"=1)
#res <- run_MCMC(startTab, data, mcmcPars, filename="test2",
#                create_post_func, mvrPars, NULL, 0.2, 
#                antigenicMap, ages, 
#                startInfHist=bestInf)

chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
#chain1 <- chain1[seq(1,nrow(chain1),by=50),]
#pairs(chain1)
plot(coda::as.mcmc(chain1))

infChain <- data.table::fread(res$history_file,data.table=FALSE)
library(plyr)
wow <- ddply(infChain, ~individual, function(x) colSums(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
wow <- reshape2::melt(wow, id.vars="individual")
wow$variable <- as.integer(wow$variable)

infHist1 <- as.data.frame(infHist)
infHist1$individual <- 1:n_indiv
colnames(infHist1) <- c(strainIsolationTimes,"individual")
infHist1 <- reshape2::melt(infHist1, id.vars="individual")
infHist1 <- infHist1[infHist1$value == 1,]

firstSamp <- infChain[infChain$sampno == 1,!(colnames(infChain) == "sampno")]
firstSamp <- reshape2::melt(firstSamp, id.vars="individual")
firstSamp <- firstSamp[firstSamp$value == 1,]

#histProfiles <- ddply(infChain, .(individual,sampno), function(x) cumsum(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
infChain <- infChain[infChain$sampno > mcmcPars["adaptive_period"],]
histProfiles <- t(apply(infChain, 1, function(x) cumsum(x[1:(ncol(infChain)-2)])))
histProfiles <- as.data.frame(cbind(histProfiles, "individual"=infChain[,c("individual")]))

histProfiles_lower <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.025))))
histProfiles_lower <- reshape2::melt(histProfiles_lower, id.vars="individual")
colnames(histProfiles_lower) <- c("individual","variable","lower")

histProfiles_upper <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.975))))
histProfiles_upper <- reshape2::melt(histProfiles_upper, id.vars="individual")
colnames(histProfiles_upper) <- c("individual","variable","upper")

histProfiles_median <- ddply(histProfiles, ~individual, function(x) apply(x, 2, function(y) quantile(y, c(0.5))))
histProfiles_median <- reshape2::melt(histProfiles_median, id.vars="individual")
colnames(histProfiles_median) <- c("individual","variable","median")

quantHist <- merge(histProfiles_lower, histProfiles_upper, by=c("individual","variable"))
quantHist <- merge(quantHist, histProfiles_median, by=c("individual","variable"))

realHistProfiles <- as.data.frame(t(apply(infHist, 1, cumsum)))
colnames(realHistProfiles) <- strainIsolationTimes
realHistProfiles$individual <- 1:n_indiv
realHistProfiles <- reshape2::melt(realHistProfiles,id.vars="individual")
  
p1 <- ggplot(quantHist[quantHist$individual %in% seq(1,10,by=1),]) + 
  geom_line(aes(x=as.integer(variable),y=median)) + 
  geom_ribbon(aes(x=as.integer(variable),ymin=lower,ymax=upper), alpha=0.2) + 
  geom_line(data=realHistProfiles[realHistProfiles$individual %in% seq(1,10,by=1),],aes(x=as.integer(variable),y=value), col="blue") +
  #scale_y_continuous(limits=c(0,10)) +
  facet_wrap(~individual, scales="free_y")

p2 <- ggplot(wow) + geom_line(aes(x=variable,y=value)) + 
    geom_vline(data=infHist1,aes(xintercept=as.integer(variable)),col="red") + 
  #geom_vline(data=firstSamp, aes(xintercept=as.integer(variable)),col="blue") +
    facet_wrap(~individual)
