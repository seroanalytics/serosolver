library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

setwd("~/Documents/Fluscape/serosolver")
#devtools::load_all()

n_indiv <- 10
fluscapeDat <- read.csv("data/fluscape_data.csv",stringsAsFactors=FALSE)
fluscapeAges <- read.csv("data/fluscape_ages.csv")
indivs <- sample(unique(fluscapeDat$individual),n_indiv)
fluscapeDat <- fluscapeDat[fluscapeDat$individual %in% indivs,]
fluscapeAges <- fluscapeAges[fluscapeAges$individual %in% indivs,]
titreDat <- fluscapeDat
ages <- fluscapeAges

titreDat$individual <- match(titreDat$individual, titreDat$individual)
ages$individual <- match(ages$individual, ages$individual)

titreDat <- titreDat[,c("individual", "samples", "virus", "titre", "run", "group")]

antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, 1)
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)

parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(1,1)
startTab <- parTab
strainIsolationTimes <- unique(fit_dat$inf_years)
samplingTimes <- unique(fluscapeDat$samples)

startInf <- setup_infection_histories_new(titreDat, strainIsolationTimes, space=5,titre_cutoff=3)

optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
f1 <- create_post_func1(optimTab,titreDat,fit_dat,NULL,infectionHistories=startInf)
startPar <- parTab$values
#startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=200))$optim$bestmem
#startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])
mcmcPars <- c("iterations"=100000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=50000,
              "save_block"=100,"thin2"=100,"histSampleProb"=0.5,"switch_sample"=2, "burnin"=0, 
              "nInfs"=4, "moveSize"=5, "histProposal"=3, "histOpt"=0)

mvrPars <- NULL

startTab$values <- startPar

#devtools::load_all()
res <- run_MCMC(startTab, titreDat, mcmcPars, filename="test1",
                create_post_func, mvrPars, PRIOR=NULL,1, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)


chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
plot(coda::as.mcmc(chain1))


infChain <- data.table::fread(res$history_file,data.table=FALSE)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
plot(coda::as.mcmc(n_hist_chain))