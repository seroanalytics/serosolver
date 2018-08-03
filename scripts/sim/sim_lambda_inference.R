library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)

## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()

## How many individuals to simulate?
n_indiv <- 69

## Which infection history proposal version to use?
describe_proposals()
histProposal <- 4

## Buckets indicates the time resolution of the analysis. Setting
## this to 1 uses annual epochs, whereas setting this to 12 gives
## monthly epochs
buckets <- 1

## The general output filename
filename <- "chains/vietnam_original_sim"

## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_lambda.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- 1
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets
parTab[parTab$names == "sigma1","values"] <- parTab[parTab$names == "sigma1","values"]
parTab[parTab$names == "sigma2","values"] <- parTab[parTab$names == "sigma2","values"]
#parTab[parTab$names == "wane","values"] <- 1/buckets
#parTab[parTab$names == "error","fixed"] <- 1
## Possible sampling times
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)

## Antigenic map for cross reactivity parameters
antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
#fit_dat <- read.csv("~/Documents/Fluscape/serosolver/data/antigenic_maps/antigenicMap_vietnam.csv")
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
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets & fit_dat$inf_years <= 2012*buckets,]
strainIsolationTimes <- unique(fit_dat$inf_years)

## Add rows for each lambda value to be inferred
tmp <- parTab[parTab$names == "lambda",]
for(i in 1:(length(strainIsolationTimes)-1)){
  parTab <- rbind(parTab, tmp)
}

## Change alpha and beta to change proposal distribution
## Setting to c(1,1) gives uniform distribution on total number of infections
#parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(2,12)

## Simulate some fake data
#strainIsolationTimes <- 1968:2015
#lambdas <- runif(length(strainIsolationTimes),0.1/buckets,0.5/buckets)
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 6*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),
                     useSIR=TRUE, pInf = NULL, useSpline=FALSE)

## If we want to use a subset of isolated strains, uncomment the line below
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets
#viruses <- seq(2000*buckets,2015*buckets,by=buckets)
#viruses <- seq(1968*buckets,2014*buckets,by=4)

titreDat <- dat[[1]]
titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]

## If we want to use pre-simulated data

#write.table(titreDat,"~/net/home/serosolver/data/sim_1000_data.csv",sep=",",row.names=FALSE)
#write.table(infHist,"~/net/home/serosolver/data/sim_1000_infHist.csv",sep=",",row.names=FALSE)
#write.table(ages,"~/net/home/serosolver/data/sim_1000_ages.csv",sep=",",row.names=FALSE)
#write.table(AR,"~/net/home/serosolver/data/sim_1000_AR.csv",sep=",",row.names=FALSE)

## Visualise simulated data
p <- plot_data(titreDat, infHist, strainIsolationTimes, 5, NULL)

## Starting infection histories based on data
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

## Housekeeping for force of infection parameters, lambda
#parTab[parTab$names == "lambda","values"] <- AR[,2]
#parTab[parTab$names == "lambda","fixed"] <- 0
#parTab[which(parTab$names == "lambda" & parTab$values <= 0),"values"] <- 0.00001


## Generate starting locations
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}

## Specify paramters controlling the MCMC procedure
mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,"adaptive_period"=50000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=1, "moveSize"=2, "histProposal"=1, "histOpt"=1,"n_par"=10)
covMat <- diag(nrow(parTab))
scale <- 0.5
w <- 1
mvrPars <- list(covMat, scale, w)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)
f <- create_post_func(parTab=startTab,data=titreDat,antigenicMap=fit_dat,PRIOR_FUNC = NULL,version=99, ageMask=ageMask)
parTab[parTab$names == "lambda","values"] <- lambdas
(f(startTab$values, infHist, 1,1,1,1))
startTab[startTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
#Rprof(tmp <- tempfile())
## Run the MCMC using the inputs generated above
#startInf <- infHist
lambdas_start <- NULL
for(year in 1:ncol(startInf)){ lambdas_start[year] <- (sum(startInf[,year])/sum(ages$DOB <= strainIsolationTimes[year]))}
startTab[startTab$names == "lambda","values"] <- lambdas_start
#parTab[parTab$names == "lambda","fixed"] <- 1
res <- run_MCMC(parTab, titreDat, mcmcPars, filename=filename,
                create_post_func, NULL,NULL,version=4, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)

#Rprof()
#summaryRprof(tmp)
#library(proftools)
#plotProfileCallGraph(readProfileData(tmp),score = "total")

#########################
## Processing outputs
#########################
## Density/trace plots
chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
plot_infection_histories(chain1, infChain, titreDat, sample(1:n_indiv, 10), fit_dat, ages,parTab,100)

y <- get_titre_predictions(chain1, infChain, titreDat, 1:n_indiv, fit_dat, ages,parTab,100, TRUE)
p1 <- ggplot(y[[3]]) + geom_histogram(aes(x=`50%`),binwidth=1) + 
  facet_wrap(~virus,scales="free_x") + 
  scale_x_continuous(limits=c(-5,5))

plot_attack_rates(infChain, titreDat, ages, 1968:2015)
wow <- infChain[infChain$i %in% 1:50,]
setkey(wow, "sampno","i")
omg <- wow[,list(V1=sum(x)), by=key(wow)]
alpha <- 100
beta <- 100
tmp <- as.data.frame(rbb(10000,48,alpha,beta))
#tmp <- data.frame(x=0:48,y=dbb(0:48,48,alpha,beta))
tmp <- data.frame(x=0:48,y=dbinom(0:48,48,prob=0.5))
ggplot(omg) + geom_density(aes(x=V1)) + geom_line(data=tmp,aes(x=x,y=y),col="red") +facet_wrap(~i)# + scale_y_continuous(limits=c(0,20))


pdf(paste0(filename, "_chain.pdf"))
plot(coda::as.mcmc(chain1))
dev.off()

## Plot inferred attack rates against true simulated attack rates
tmp <- summary(as.mcmc(chain1[,2:(ncol(chain1)-1)]))
tmp <- as.data.frame(tmp[[2]])
tmp$names <- rownames(tmp)
lambda_names <- c("lambda",paste0("lambda.",1:(nrow(AR)-1)))
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
 # scale_x_continuous(breaks=seq(1968,2015,by=1)) +
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
indivs <- sample(n_indiv, 10)
sampd <- sample(n_indiv,20)
infChain <- data.table::fread(res$history_file)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
n_strain <- max(infChain$j)
data.table::setkey(infChain, "i","sampno")
n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
inf_chain_p <- ggplot(n_inf_chain[n_inf_chain$i %in% sampd,]) + geom_line(aes(x=sampno,y=V1)) + facet_wrap(~i)
## Generate cumulative infection history plots for a random subset of individuals
## based on data and posterior
#mcmcPars["adaptive_period"]+mcmcPars["burnin"]

ps <- generate_cumulative_inf_plots(res$history_file, mcmcPars["adaptive_period"], 
                                           sampd, infHist, startInf,strainIsolationTimes, 100,ages,numberCol=4)
plot_infection_histories(chain1, infChain, titreDat, sample(1:100, 10), fit_dat, ages,parTab,100)
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



