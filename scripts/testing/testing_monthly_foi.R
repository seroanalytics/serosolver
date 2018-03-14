library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
## Set working directory and load code
setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
n_indiv <-200
buckets <- 12
nYears <- length(1968:2015)
error <- 0.05
foi <- runif(nYears, 0.1,0.4)
x <- seq(0,buckets-1,by=1)/buckets
knots <- c(0.33,0.66)
degree <- 2
nKnots <- length(knots) + degree + 1
theta <- runif(nKnots)

filename <- "chains/test_foi"
## Read in parameter table to simulate from and change waning rate if necessary
parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
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
## All possible circulation times
fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets,]
strainIsolationTimes <- unique(fit_dat$inf_years)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 10*buckets,75*buckets,
                     simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),knots=knots,theta=theta,
                     useSIR=FALSE,useSpline=TRUE)
viruses <- c(1968, 1969, 1972, 1975, 1977, 1979, 1982, 1985, 1987, 
             1989, 1992, 1995, 1998, 2000, 2002, 2004, 2007, 2009, 
             2010, 2012, 2014)*buckets
titreDat <- dat[[1]]
titreDat <- titreDat[titreDat$virus %in% viruses,]
infectionHistories <- infHist <- dat[[2]]
ages <- dat[[3]]
AR <- dat[[4]]
lambdas <- dat[[5]]
p <- plot_data(titreDat, infHist, strainIsolationTimes, 5, NULL)
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

parTab$identity <- 1
knotsTab <- data.frame(names=c("knot1","knot2"),
                       values=knots,
                       fixed=0,steps=0.01,lower_bound=0,upper_bound=1,lower_start=0,upper_start=1,
                       identity=4)
thetaTab <- data.frame(names=paste0("theta",1:length(theta)),
                       values=theta,
                       fixed=0,steps=0.01,lower_bound=0,upper_bound=1,lower_start=0,upper_start=1,
                       identity=3)
foiTab <- data.frame(names=paste0("foi",1:length(foi)),
                       values=foi,
                       fixed=0,steps=0.01,lower_bound=0,upper_bound=1,lower_start=0,upper_start=1,
                       identity=2)
parTab <- rbind(parTab,foiTab,thetaTab,knotsTab)

f <- create_post_func(parTab,titreDat,fit_dat,NULL,7,ageMask)
sum(f(parTab$values, infHist))

## Generate starting locations
startTab <- parTab
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}
mcmcPars <- c("iterations"=1000000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=5000,"thin"=10,"adaptive_period"=500000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=5, "burnin"=0, 
              "nInfs"=4, "moveSize"=2, "histProposal"=3, "histOpt"=1)
res <- run_MCMC(startTab, titreDat, mcmcPars, filename=filename,
                create_post_func, NULL, NULL,version=7, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)

chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
pdf(paste0(filename, "_chain.pdf"))
plot(coda::as.mcmc(chain1))
dev.off()

tmp <- summary(as.mcmc(chain1))
tmp <- as.data.frame(tmp[[2]])
tmp$names <- rownames(tmp)
lambda_names <- paste0("foi",1:48)
tmp <- tmp[tmp$names %in% lambda_names,]
tmpTab <- parTab[parTab$identity == 2,]
tmpTab$values <- zikaProj::sum_buckets(dat[[5]],rep(12,48))
tmpTab$names <- lambda_names
#tmp$names <- strainIsolationTimes
#tmpTab$names <- strainIsolationTimes
tmp$names <- 1968:2015
tmpTab$names <- 1968:2015
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
