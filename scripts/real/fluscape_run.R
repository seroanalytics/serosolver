library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

## How many individuals to fit to?
n_indiv <-100
FLUSCAPE <- TRUE
LAMBDA <- TRUE

## Which infection history proposal version to use?
## ********Note*********
## that this should be set to 1 if using the explicit FOI term
describe_proposals()
histProposal <- 6



## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1

## The general output filename
filename <- "testing/fluscape"

## Read in parameter table to simulate from and change waning rate if necessary
## Use "parTab.csv" if not using explicit phi term, "parTab_phis.csv" otherwise
if(LAMBDA){
    parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_phi.csv",stringsAsFactors=FALSE)
} else {
    parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
}

## Possible sampling times
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)

## Antigenic map for cross reactivity parameters
fit_dat <- read.csv("~/Documents/Fluscape/serosolver/data/antigenic_maps/antigenicMap_vietnam.csv")

## Read in Fluscape or Vietnam data
if(FLUSCAPE){
    fluscapeDat <- read.csv("data/real/fluscape_data.csv",stringsAsFactors=FALSE)
    fluscapeAges <- read.csv("data/real/fluscape_ages.csv")
    ## Remove individuals with NA for DOB
    na_indiv <- fluscapeAges[which(is.na(fluscapeAges$DOB)),"individual"]
    fluscapeDat <- fluscapeDat[-na_indiv,]
    fluscapeAges <- fluscapeAges[-na_indiv,]
    ## Take random subset of individuals
    indivs <- sample(unique(fluscapeDat$individual),n_indiv)
    indivs <- indivs[order(indivs)]
    #indivs <- unique(fluscapeAges$individual)
    titreDat <- fluscapeDat[fluscapeDat$individual %in% indivs,]
    ages <- fluscapeAges[fluscapeAges$individual %in% indivs,]
    titreDat$individual <- match(titreDat$individual, indivs)
    ages$individual <- match(ages$individual, indivs)
    titreDat <- titreDat[,c("individual", "samples", "virus", "titre", "run", "group")]
} else {
    titreDat <- read.csv("data/real/vietnam_data.csv")
    ages <- data.frame(individual=1:length(unique(titreDat$individual)),DOB=1940)
    fit_dat <- fit_dat[fit_dat$inf_years <= 2012,]
}

## All possible circulation times
strainIsolationTimes <- unique(fit_dat$inf_years)

##############
## BETA BINOMIAL PROPOSALS
## NOTE - this is not relevant if using history proposal 1
##############
## Change alpha and beta to change proposal distribution
## Setting to c(1,1) gives uniform distribution on total number of infections
#parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)


## Starting infection histories based on data
## Can use starting history using the algorithm in the bioarxiv paper or a new, sparser one
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
#startInf <- matrix(sample(c(0,0,0,0,0,0,0,0,1),n_indiv*length(strainIsolationTimes),replace=TRUE),nrow=n_indiv)
#startInf <- setup_infection_histories_OLD(titreDat, unique(fit_dat$inf_years), rep(1,n_indiv), sample_prob=0.2, titre_cutoff=3)

if(LAMBDA){
    if(!FLUSCAPE){
        ages1 <- read.csv("~/Documents/Fluscape/serosolver/data/real/vietnam_ages.csv")
    } else {
        ages1 <- ages
    }
    
    version <- 4
    n_alive <- sapply(strainIsolationTimes, function(x) length(ages[ages$DOB <= x,])/69)
    tmp <- parTab[parTab$names == "phi",]
    for(i in 1:(length(strainIsolationTimes)-1)){
        parTab <- rbind(parTab, tmp)
    }
    parTab[parTab$names == "phi","upper_bound"] <- n_alive
    parTab[parTab$names == "phi","upper_start"] <- n_alive
} else {
    version <- 1
}
parTab <- parTab[parTab$names != "phi",]
## Generate starting parameters for the MCMC chain
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

## Multivariate proposals or univariate? Use univariate for now
covMat <- diag(nrow(parTab))
scale <- 0.005
w <- 1
mvrPars <- list(covMat, scale, w)
mvrPars <- NULL

#########################
## RUN MCMC
#########################
startTab[startTab$names == "wane","fixed"] <- 0
startTab[startTab$names == "wane","values"] <- 1
startTab[startTab$names == "wane","upper_bound"] <- 1

## MCMC run time
mcmcPars <- c("iterations"=70000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,"adaptive_period"=20000,
              "save_block"=1000,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=floor(ncol(startInf)/2), "moveSize"=5, "histProposal"=6, "histOpt"=0,"n_par"=10, "swapPropn"=0.5)


res <- run_MCMC(startTab, titreDat[titreDat$run == 1,], mcmcPars, filename=filename,
                create_post_func, mvrPars, PRIOR=NULL,version=1, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)
#########################


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


n_strain <- max(infChain$j)
data.table::setkey(infChain, "i","sampno")
n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
n_inf_chain$chain <- i
ess_inf <- sapply(unique(n_inf_chain$i), function(x) effectiveSize(n_inf_chain[n_inf_chain$i == x,"V1"]))
## Order by increasing ESS
ess_worst <- order(ess_inf)
summary(ess_inf)
## Which individuals to look at
indivs <- c(ess_worst[1:12])
ggplot(n_inf_chain[n_inf_chain$i %in% indivs, ]) + 
  geom_line(aes(x=sampno, y=V1,col=as.factor(chain))) + 
  theme_bw() +
  ylab("Total number of infections") +
  xlab("MCMC iteration") +
  scale_y_continuous(limits=c(0,45),expand=c(0,0))+
  facet_wrap(~i)

ggplot(n_inf_chain[n_inf_chain$i %in% indivs, ]) + 
  geom_density(aes(x=V1,col=as.factor(chain))) + 
  theme_bw() +
  facet_wrap(~i)


indivHist <- plyr::ddply(n_inf_chain,.(i),function(x) quantile(x$V1, c(0,0.025,0.5,0.975,1)))
colnames(indivHist) <- c("individual","min","lower","median","upper","max")
indivHist <- indivHist[order(indivHist$median),]
indivHist$individual <- 1:nrow(indivHist)
p <- ggplot(indivHist) + 
  geom_pointrange(aes(x=individual,y=median,ymin=lower,ymax=upper),
                  size=0.1,shape=21,fatten=0.1) +
  scale_y_continuous(limits=c(0,45),expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,1100,by=100)) +
  xlab("Individual") +
  ylab("Estimated number of infections") +
  theme_bw()




if(FLUSCAPE){
    n_alive <- NULL
} else {
    n_alive <- sapply(strainIsolationTimes, function(x) length(ages1[ages1$DOB <= x,]))
}
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
#generate_all_plots(getwd(), mcmcPars["adaptive_period"] + mcmcPars["burnin"], res$chain_file, res$history_file,
#                   titreDat, fit_dat, parTab, ages, nIndiv=10,nSamp=100,
#                   filename=filename)

