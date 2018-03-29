library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

## How many individuals to fit to?
#n_indiv <-10

## Which infection history proposal version to use?
describe_proposals()
histProposal <- 1

## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1

## The general output filename
filename <- "vietnam"

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets

## Possible sampling times
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)

## Antigenic map for cross reactivity parameters
fit_dat <- read.csv("~/Documents/Fluscape/serosolver/data/antigenicMap_AK.csv")

## Read in Fluscape data
fluscapeDat <- read.csv("data/fluscape_data.csv",stringsAsFactors=FALSE)
fluscapeAges <- read.csv("data/fluscape_ages.csv")

titreDat <- read.csv("data/vietnam.csv")
ages <- data.frame(individual=1:length(unique(titreDat$individual)),DOB=1940)
fit_dat <- fit_dat[fit_dat$inf_years <= 2012,]

## Remove individuals with NA for DOB
na_indiv <- fluscapeAges[which(is.na(fluscapeAges$DOB)),"individual"]
fluscapeDat <- fluscapeDat[-na_indiv,]
fluscapeAges <- fluscapeAges[-na_indiv,]

## Take random subset of individuals
#indivs <- sample(unique(fluscapeDat$individual),n_indiv)
#indivs <- indivs[order(indivs)]
indivs <- unique(fluscapeAges$individual)
titreDat <- fluscapeDat[fluscapeDat$individual %in% indivs,]
ages <- fluscapeAges[fluscapeAges$individual %in% indivs,]
titreDat$individual <- match(titreDat$individual, indivs)
ages$individual <- match(ages$individual, indivs)
titreDat <- titreDat[,c("individual", "samples", "virus", "titre", "run", "group")]

## All possible circulation times
strainIsolationTimes <- unique(fit_dat$inf_years)

## Change alpha and beta to change proposal distribution
## Setting to c(1,1) gives uniform distribution on total number of infections
#parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(0.75,4.5)

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=3)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

## Generate starting locations for the parameter vector, theta.
## Generate this by optimising theta based on the chosen starting infection histories
startTab <- parTab
optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
f1 <- create_post_func1(optimTab,titreDat,fit_dat,NULL,infectionHistories=startInf)
startPar <- parTab$values
startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=100))$optim$bestmem
startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])
startTab$values <- startPar

## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=1000,"thin"=1,"adaptive_period"=50000,
              "save_block"=100,"thin2"=10,"histSampleProb"=0.5,"switch_sample"=2, "burnin"=0, 
              "nInfs"=1, "moveSize"=10, "histProposal"=4, "histOpt"=0)
startTab[startTab$names == "wane","values"] <- 1
titreDat1 <- titreDat[titreDat$run == 1,]
startTab[startTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
system.time(
  
## Run the MCMC using the inputs generated above
res <- run_MCMC(startTab, titreDat, mcmcPars, filename=filename,
                create_post_func, NULL, PRIOR,version=1, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)
)

#########################
## Processing outputs
#########################
## Density/trace plots
chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
pdf(paste0(filename, "_chain.pdf"))
plot(coda::as.mcmc(chain1))
dev.off()

## Plot inferred attack rates
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
xs <- min(strainIsolationTimes):max(strainIsolationTimes)
colnames(AR) <- c("year","AR")
arP <- plot_attack_rates(infChain, titreDat,ages,xs, n_alive) + scale_y_continuous(limits=c(0,1), expand=c(0,0))


## Density/trace plots on total number of infections
# The ddply call takes a while, so comment back in if you are willing to wait
#n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
#n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
#n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
#pdf(paste0(filename,"_infChain.pdf"))
#plot(coda::as.mcmc(n_hist_chain))
#dev.off()

## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
infHist_p <- generate_cumulative_inf_plots(res$history_file, mcmcPars["adaptive_period"]+mcmcPars["burnin"], 10, NULL)
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

