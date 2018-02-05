library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

## How many individuals to simulate?
n_indiv <-10

## Which infection history proposal version to use?
describe_proposals()
histProposal <- 2

## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1

## The general output filename
filename <- "test_monthly"

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets

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

## Change alpha and beta to change proposal distribution
## Setting to c(1,1) gives uniform distribution on total number of infections
#parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)

## Simulate some fake data
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 10*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=FALSE)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)

titreDat <- dat[[1]]
##titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

## Code to save simulated data
##write.table(titreDat,"~/net/home/serosolver/data/sim_200_dat.csv",row.names=FALSE,sep=",")
##write.table(infHist,"~/net/home/serosolver/data/sim_200_infHist.csv",row.names=FALSE,sep=",")
##write.table(ages,"~/net/home/serosolver/data/sim_200_ages.csv",row.names=FALSE,sep=",")
##write.table(AR,"~/net/home/serosolver/data/sim_200_AR.csv",row.names=FALSE,sep=",")

## If we want to use pre-simulated data
#titreDat <- read.csv("data/sim_200_titres.csv",stringsAsFactors=FALSE)
#ages <- read.csv("data/sim_200_ages.csv",stringsAsFactors=FALSE)
#infectionHistories <- infHist <- read.csv("data/sim_200_infHist.csv",stringsAsFactors=FALSE)
#AR <- read.csv("data/sim_200_AR.csv",stringsAsFactors=FALSE)

## Visualise simulated data
p <- plot_data(titreDat, infHist, strainIsolationTimes, 5, NULL)

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

## Generate starting locations for the parameter vector, theta.
## Generate this by optimising theta based on the chosen starting infection histories
startTab <- parTab
optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
f1 <- create_post_func1(optimTab,titreDat,fit_dat,NULL,infectionHistories=startInf)
startPar <- parTab$values
startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=10))$optim$bestmem
startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])
startTab$values <- startPar

## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=20000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=1000,"thin"=1,"adaptive_period"=10000,
              "save_block"=100,"thin2"=1,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=4, "moveSize"=5, "histProposal"=histProposal, "histOpt"=0)

## Run the MCMC using the inputs generated above
res <- run_MCMC(startTab, titreDat, mcmcPars, filename=filename,
                create_post_func, NULL, PRIOR,version=1, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)

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
infChain <- data.table::fread(res$history_file,data.table=FALSE)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
xs <- min(strainIsolationTimes):max(strainIsolationTimes)
colnames(AR) <- c("year","AR")
arP <- plot_attack_rates(infChain, titreDat,ages,xs) + geom_point(data=AR,aes(x=year,y=AR), col="green")


## Density/trace plots on total number of infections
#n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
#n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
#n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
#pdf(paste0(filename,"_infChain.pdf"))
#plot(coda::as.mcmc(n_hist_chain))
#dev.off()

## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
infHist_p <- generate_cumulative_inf_plots(res$history_file, mcmcPars["adaptive_period"]+mcmcPars["burnin"], 10, infHist)
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

