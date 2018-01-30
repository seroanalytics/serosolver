library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

setwd("~/Documents/serosolver_own/serosolver")
devtools::load_all()

filename <- "chains/output_multiprob"
n_indiv <-50

fluscapeDat <- read.csv("data/fluscape_data.csv",stringsAsFactors=FALSE)
fluscapeAges <- read.csv("data/fluscape_ages.csv")
na_indiv <- fluscapeAges[which(is.na(fluscapeAges$DOB)),"individual"]
fluscapeDat <- fluscapeDat[-na_indiv,]
fluscapeAges <- fluscapeAges[-na_indiv,]
#indivs <- 1:10
indivs <- sample(unique(fluscapeDat$individual),n_indiv)
indivs <- indivs[order(indivs)]
titreDat <- fluscapeDat[fluscapeDat$individual %in% indivs,]
ages <- fluscapeAges[fluscapeAges$individual %in% indivs,]

titreDat$individual <- match(titreDat$individual, indivs)
ages$individual <- match(ages$individual, indivs)
>>>>>>> cfa9d030ac0ca384f5a6658a792ace3715ab14af

titreDat <- titreDat[,c("individual", "samples", "virus", "titre", "run", "group")]

#antigenicMap <- read.csv("~/Documents/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
#fit_dat <- generate_antigenic_map(antigenicMap, 1)
fit_dat <- read.csv("~/Documents/serosolver_own/serosolver/data/antigenicMap_AK.csv")
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)

parTab <- read.csv("~/Documents/serosolver_own/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(2,12)
startTab <- parTab
strainIsolationTimes <- unique(fit_dat$inf_years)
samplingTimes <- unique(fluscapeDat$samples)

startInf <- setup_infection_histories_new(titreDat, ages, strainIsolationTimes, space=4,titre_cutoff=3)

optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
f1 <- create_post_func1(optimTab,titreDat,fit_dat,NULL,infectionHistories=startInf)
startPar <- parTab$values
startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=20))$optim$bestmem
startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])
mcmcPars <- c("iterations"=100000,"popt"=0.234,"popt_hist"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=50000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
              "nInfs"=3, "moveSize"=2, "histProposal"=3, "histOpt"=1)

covMat <- diag(nrow(parTab))
scale <- 0.01
w <- 0.95
mvrPars <- list(covMat, scale, w)
mvrPars <- NULL

#startPar <- c(1.2,3.3,0.03,0.2,0.1,0.03,8,1,1,4.5)
startTab$values <- startPar
startTab[startTab$names == "wane","values"] <- 0.2
#devtools::load_all()
res <- run_MCMC(startTab, titreDat, mcmcPars, filename=filename,
                create_post_func, mvrPars, PRIOR=NULL,1, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)


chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
pdf("testing_real.pdf")
plot(coda::as.mcmc(chain1))
dev.off()

antigenicMap <- read.csv("~/Documents/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
generate_all_plots(getwd(), mcmcPars["adaptive_period"] + mcmcPars["burnin"], res$chain_file, res$history_file,
                            titreDat, fit_dat, parTab, ages, nIndiv=10,nSamp=100,
                               filename=filename)

infChain <- data.table::fread(res$history_file,data.table=FALSE)
#infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
pdf("testing_real_hist.pdf")
plot(coda::as.mcmc(n_hist_chain))
dev.off()






infChain <- data.table::fread(res$history_file,data.table=FALSE)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
wow <- plyr::ddply(infChain, ~individual, function(x) colSums(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
wow <- reshape2::melt(wow, id.vars="individual")
wow$variable <- as.integer(wow$variable)
firstSamp <- infChain[infChain$sampno == 1,!(colnames(infChain) == "sampno")]
firstSamp <- reshape2::melt(firstSamp, id.vars="individual")
firstSamp <- firstSamp[firstSamp$value == 1,]

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

startHistProfiles <- as.data.frame(t(apply(startInf, 1, cumsum)))
colnames(startHistProfiles) <- strainIsolationTimes
startHistProfiles$individual <- 1:n_indiv
startHistProfiles <- startHistProfiles[startHistProfiles$individual %in% 1:10,]
startHistProfiles <- reshape2::melt(startHistProfiles,id.vars="individual")

p1 <- ggplot(quantHist[quantHist$individual %in% seq(1,10,by=1),]) + 
  geom_line(aes(x=as.integer(variable),y=median)) + 
  geom_ribbon(aes(x=as.integer(variable),ymin=lower,ymax=upper), alpha=0.2) + 
  geom_line(data=startHistProfiles[startHistProfiles$individual %in% seq(1,10,by=1),],aes(x=as.integer(variable),y=value), col="red") +
  facet_wrap(~individual, scales="free_y")



svg(paste0(filename,"_infectionhistories.svg"))
plot(p1)
dev.off()

