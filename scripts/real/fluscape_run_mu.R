library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

## How many individuals to fit to?
n_indiv <-1000

## Which infection history proposal version to use?
describe_proposals()
histProposal <- 1

## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1
mu_bucket <- 1

## The general output filename
filename <- "fluscape_3"

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_mus.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets

## Possible sampling times
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)

## Antigenic map for cross reactivity parameters
#fit_dat <- read.csv("~/Documents/Fluscape/serosolver/data/antigenicMap_AK.csv")
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)

## Read in Fluscape data
fluscapeDat <- read.csv("data/fluscape_data.csv",stringsAsFactors=FALSE)
fluscapeAges <- read.csv("data/fluscape_ages.csv")

## Remove individuals with NA for DOB
na_indiv <- fluscapeAges[which(is.na(fluscapeAges$DOB)),"individual"]
fluscapeDat <- fluscapeDat[!(fluscapeDat$individual %in% na_indiv),]
fluscapeAges <- fluscapeAges[!(fluscapeAges$individual %in% na_indiv),]

## Take random subset of individuals
#indivs <- sample(unique(fluscapeDat$individual),n_indiv)
indivs<- intersect(unique(fluscapeDat$individual),unique(fluscapeAges$individual))
indivs <- indivs[order(indivs)]
#indivs <- unique(fluscapeAges$individual)

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
n_mus <- length(strainIsolationTimes)/mu_bucket
mu_indices <- rep(1:n_mus, each=mu_bucket) - 1

## Add rows for each lambda value to be inferred
tmp <- parTab[parTab$names == "lambda",]
for(i in 1:(length(strainIsolationTimes)-1)){
  parTab <- rbind(parTab, tmp)
}

## Generate value for each mu
mu_mean <- parTab[parTab$names == "mu_mean","values"]
mu_sd <- parTab[parTab$names == "mu_sd","values"]
l_mean <- log(mu_mean) - (mu_sd^2)/2
mus <- rlnorm(n_mus, mean=l_mean, sd=mu_sd)
mu_tab <- data.frame(names="mu",values=mus,fixed=0,steps=0.1,lower_bound=0,upper_bound=8,lower_start=0.5,
                     upper_start=5,identity=3,block=1)
parTab <- rbind(mu_tab, parTab)

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=3)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

## Generate starting locations for the parameter vector, theta.
## Generate this by optimising theta based on the chosen starting infection histories
## Generate starting locations
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

#startTab <- parTab
#optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
#f1 <- create_post_func1(optimTab,titreDat,fit_dat,NULL,infectionHistories=startInf)
#startPar <- parTab$values
#startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=100))$optim$bestmem
#startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])
#startTab$values <- startPar

## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=2000000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=5000,"thin"=10,"adaptive_period"=500000,
              "save_block"=1000,"thin2"=500,"histSampleProb"=1,"switch_sample"=2, "burnin"=50000, 
              "nInfs"=4, "moveSize"=5, "histProposal"=1, "histOpt"=1,"n_par"=10)

system.time(
## Run the MCMC using the inputs generated above
res <- run_MCMC(startTab, titreDat, mcmcPars, filename=filename,
                create_post_func_mu, NULL, PRIOR,version=4, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf, mu_indices=mu_indices)
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
arP <- plot_attack_rates(infChain, titreDat,ages,xs)
sampd <- sample(1:n_indiv, 10)
plot_infection_histories(chain1, infChain, titreDat, sampd, 
                         fit_dat, ages,parTab,100, mu_indices-1)

## Density/trace plots on total number of infections
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
n_strain <- max(infChain$j)
data.table::setkey(infChain, "i","sampno")
n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
inf_chain_p <- ggplot(n_inf_chain[n_inf_chain$i %in% sampd,]) + geom_line(aes(x=sampno,y=V1)) + facet_wrap(~i)

## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
ps <- generate_cumulative_inf_plots(res$history_file, mcmcPars["adaptive_period"], 
                                    sampd, startInf, startInf,strainIsolationTimes, 100,ages)
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

