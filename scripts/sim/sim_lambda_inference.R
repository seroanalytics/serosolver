library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

## How many individuals to simulate?
n_indiv <-1000

## Which infection history proposal version to use?
describe_proposals()
histProposal <- 1

## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 12

## The general output filename
filename <- "chains/test_lambda"

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_lambda.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets
parTab[parTab$names == "wane","values"] <- 1
## Possible sampling times
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)

## Antigenic map for cross reactivity parameters
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)

## Rename circulation years based on isolation time
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets
antigenicMap$Strain <- virus_key[antigenicMap$Strain]

## For visualisation - the antigenic summary path
p1 <- ggplot(antigenicMap) + 
  geom_line(data=fit_dat,aes(x=x_coord,y=y_coord), col="red") +
  geom_point(data=antigenicMap,aes(x=X,y=Y)) + 
  geom_label(data=antigenicMap,aes(x=X+4,y=Y+0.25,label=Strain)) +
  theme_bw()

## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets,]
strainIsolationTimes <- unique(fit_dat$inf_years)

## Add rows for each lambda value to be inferred
tmp <- parTab[parTab$names == "lambda",]
for(i in 1:(length(strainIsolationTimes)-1)){
  parTab <- rbind(parTab, tmp)
}

## Change alpha and beta to change proposal distribution
## Setting to c(1,1) gives uniform distribution on total number of infections
#parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)

## Simulate some fake data
#strainIsolationTimes <- 1968:2015
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 10*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),
                     useSIR=TRUE)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets

titreDat <- dat[[1]]
titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

## If we want to use pre-simulated data
#titreDat <- read.csv("data/sim_1000_titres.csv",stringsAsFactors=FALSE)
#ages <- read.csv("data/sim_1000_ages.csv",stringsAsFactors=FALSE)
#infectionHistories <- infHist <- read.csv("data/sim_1000_infHist.csv",stringsAsFactors=FALSE)
#AR <- read.csv("data/sim_1000_AR.csv",stringsAsFactors=FALSE)

#write.table(titreDat,"data/sim_10_titres.csv",sep=",",row.names=FALSE)
#write.table(infHist,"data/sim_10_infHist.csv",sep=",",row.names=FALSE)
#write.table(ages,"data/sim_10_ages.csv",sep=",",row.names=FALSE)
#write.table(AR,"data/sim_10_AR.csv",sep=",",row.names=FALSE)

## Visualise simulated data
p <- plot_data(titreDat, infHist, strainIsolationTimes, 5, NULL)

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

## Housekeeping for force of infection parameters, lambda
parTab[parTab$names == "lambda","values"] <- AR[,2]
parTab[parTab$names == "lambda","fixed"] <- 0
parTab[which(parTab$names == "lambda" & parTab$values <= 0),"values"] <- 0.00001


## Generate starting locations
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=500000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=100000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=10, "burnin"=0, 
              "nInfs"=4, "moveSize"=2, "histProposal"=histProposal, "histOpt"=1)
covMat <- diag(nrow(parTab))
scale <- 0.01
w <- 0.9
mvrPars <- list(covMat, scale, w)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)
f <- create_post_func(parTab=startTab,data=titreDat,antigenicMap=fit_dat,PRIOR_FUNC = NULL,version=4, ageMask=ageMask)
f(startTab$values, startInf)
Rprof(tmp <- tempfile())
## Run the MCMC using the inputs generated above
res <- run_MCMC(startTab, titreDat, mcmcPars, filename=filename,
                create_post_func, NULL, NULL,version=4, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)

Rprof()
summaryRprof(tmp)
library(proftools)
plotProfileCallGraph(readProfileData(tmp),score = "total")

#########################
## Processing outputs
#########################
## Density/trace plots
chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
pdf(paste0(filename, "_chain.pdf"))
plot(coda::as.mcmc(chain1))
dev.off()

## Plot inferred attack rates against true simulated attack rates
tmp <- summary(as.mcmc(chain1))
tmp <- as.data.frame(tmp[[2]])
tmp$names <- rownames(tmp)
lambda_names <- c("lambda",paste0("lambda.",1:47))
tmp <- tmp[tmp$names %in% lambda_names,]
tmpTab <- parTab[parTab$names == "lambda",]
tmpTab$values <- AR[,2]
tmpTab$names <- lambda_names
tmp$names <- strainIsolationTimes
tmpTab$names <- strainIsolationTimes
AR_recovery <- ggplot() +
  geom_pointrange(data=tmp,aes(x=names,y=`50%`,ymin=`2.5%`,ymax=`97.5%`)) +
  geom_point(data=tmpTab,aes(x=names,y=values),col="red") +
  ylab("Attack rate, lambda") +
  xlab("Year") +
  scale_y_continuous(limits=c(0,1)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1))
svg(paste0(filename,"actual_AR.svg"))
plot(AR_recovery)
dev.off()

## Density/trace plots on total number of infections
#n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
#n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
#n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
#pdf(paste0(filename,"_infChain.pdf"))
#plot(coda::as.mcmc(n_hist_chain))
#dev.off()

## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
infHist_p <- generate_cumulative_inf_plots(res$history_file, mcmcPars["adaptive_period"]+mcmcPars["burnin"], 
                                           10, infHist, startInf,strainIsolationTimes)
svg(paste0(filename, "cumulative_infHist.svg"))
plot(infHist_p)
dev.off()

## See help file for this function for details about plot outputs
## 1. Correlation plots
## 2. Autocorrelation plots
## 3. Posterior density plots
## 4. Posterior trace plots
## 5. Model fits over titre data
## 6. Attack rate plots, as above
## 7. Total number of infections for all individuals
generate_all_plots(getwd(), mcmcPars["adaptive_period"] + mcmcPars["burnin"], res$chain_file, res$history_file,
                   titreDat, fit_dat, parTab, ages, nIndiv=10,nSamp=100,
                   filename=filename)

