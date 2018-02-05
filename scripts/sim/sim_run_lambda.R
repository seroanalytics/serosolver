library(ggplot2)
library(coda)
library(plyr)
library(reshape2)

setwd("~/Documents/Fluscape/serosolver")
devtools::load_all()
n_indiv <- 200

version <- 4
histProposal <- 3
PRIOR <- infHistPrior
buckets <- 1

antigenicMap <- read.csv("~/Documents/Fluscape/fluscape/trunk/data/Fonville2014AxMapPositionsApprox.csv",stringsAsFactors=FALSE)
fit_dat <- generate_antigenic_map(antigenicMap, buckets)
virus_key <- c("HK68"=1968, "EN72"=1972, "VI75"=1975, "TX77"=1977, "BK79"=1979, "SI87"=1987, "BE89"=1989, "BJ89"=1989,
               "BE92"=1992, "WU95"=1995, "SY97"=1997, "FU02"=2002, "CA04"=2004, "WI05"=2005, "PE06"=2006)*buckets
antigenicMap$Strain <- virus_key[antigenicMap$Strain]

p1 <- ggplot(antigenicMap) + 
  geom_line(data=fit_dat,aes(x=x_coord,y=y_coord), col="red") +
  geom_point(data=antigenicMap,aes(x=X,y=Y)) + 
  geom_label(data=antigenicMap,aes(x=X+4,y=Y+0.25,label=Strain)) +
  theme_bw()


fit_dat <- fit_dat[fit_dat$inf_years >= 1968*buckets,]
strainIsolationTimes <- unique(fit_dat$inf_years)
samplingTimes <- seq(2010*buckets, 2015*buckets, by=1)


parTab <- read.csv("~/Documents/Fluscape/serosolver/inputs/parTab_lambda.csv",stringsAsFactors=FALSE)
parTab[parTab$names == "wane","values"] <- parTab[parTab$names == "wane","values"]/buckets
tmp <- parTab[parTab$names == "lambda",]
for(i in 1:(length(strainIsolationTimes)-1)){
  parTab <- rbind(parTab, tmp)
}
parTab[parTab$names %in% c("alpha","beta"),"values"] <- find_a_b(length(strainIsolationTimes),7,50)
parTab[parTab$names %in% c("alpha","beta"),"values"] <- c(100,100)
parTab[parTab$names == "wane","values"] <- 0.5


titreDat <- read.csv("data/sim_200_lambda_titres.csv",stringsAsFactors=FALSE)
ages <- read.csv("data/sim_200_lambda_ages.csv",stringsAsFactors=FALSE)
infectionHistories <- infHist <- read.csv("data/sim_200_lambda_infHist.csv",stringsAsFactors=FALSE)
AR <- read.csv("data/sim_200_lambda_AR.csv",stringsAsFactors=FALSE)

#dat <- simulate_data(parTab, 1, n_indiv, buckets,strainIsolationTimes,
#                     samplingTimes, 2, antigenicMap=fit_dat, 0, 0, 10*buckets,75*buckets,
 #                    simInfPars=c("mean"=0.15,"sd"=0.5,"bigMean"=0.5,"logSD"=1),useSIR=FALSE)
#titreDat <- dat[[1]]
#viruses <- c(1968L, 1969L, 1972L, 1975L, 1977L, 1979L, 1982L, 1985L, 1987L, 
#  1989L, 1992L, 1995L, 1998L, 2000L, 2002L, 2004L, 2007L, 2009L, 
#  2010L, 2012L, 2014L)
#titreDat <- titreDat[titreDat$virus %in% viruses,]
#infectionHistories <- infHist <- dat[[2]]
#ages <- dat[[3]]
#AR <- dat[[4]]
p <- plot_data(titreDat, infHist, strainIsolationTimes, 5, NULL)

parTab[parTab$names == "lambda","values"] <- AR[,2]
parTab[parTab$names == "lambda","fixed"] <- 0
parTab[which(parTab$names == "lambda" & parTab$values <= 0),"values"] <- 0.00001
startTab <- parTab
startInf <- setup_infection_histories_new(titreDat, ages, unique(fit_dat$inf_years), space=5,titre_cutoff=2)
ageMask <- create_age_mask(ages, strainIsolationTimes,n_indiv)

#optimTab <- startTab[!(startTab$names %in% c("alpha","beta")),]
#f1 <- create_post_func1(optimTab,titreDat,fit_dat,NULL,infectionHistories=startInf)
#startPar <- parTab$values
#startPar <- DEoptim::DEoptim(f1, lower=optimTab$lower_bound, upper=optimTab$upper_bound,control=list(itermax=10))$optim$bestmem
#startPar <- c(startPar, startTab[(startTab$names %in% c("alpha","beta")),"values"])
#startTab[!(startTab$names %in% c("alpha","beta")),"values"] <- startPar
#startPar[6] <- 0.1
for(i in 1:nrow(startTab)){
  if(startTab[i,"fixed"] == 0){
    startTab[i,"values"] <- runif(1,startTab[i,"lower_start"], 
                                  startTab[i,"upper_start"])
  }
}
mcmcPars <- c("iterations"=2000000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=10,"adaptive_period"=500000,
              "save_block"=100,"thin2"=100,"histSampleProb"=1,"switch_sample"=10, "burnin"=0, 
              "nInfs"=4, "moveSize"=5, "histProposal"=histProposal, "histOpt"=1)

## For univariate proposals
covMat <- diag(nrow(parTab))
scale <- 0.01
w <- 0.95
mvrPars <- list(covMat, scale, w)

mvrPars <- NULL
#devtools::load_all()
f <- create_post_func(parTab,titreDat,fit_dat, NULL,4,ageMask)
#mvrPars <- NULL
#f(parTab$values, infHist)
res <- run_MCMC(startTab, titreDat, mcmcPars, filename="lambda_sim_200",
                create_post_func, mvrPars, PRIOR,version=version, 0.2, 
                fit_dat, ages=ages, 
                startInfHist=startInf)


chain1 <- read.csv(res$chain_file)
chain1 <- chain1[chain1$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
pdf("lammbda_sim_200_chain.pdf")
plot(coda::as.mcmc(chain1))
dev.off()

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
    theme_bw() +
    theme(axis.text.x=element_text(angle=45,hjust=1))
svg("lambda_sim_200_AR.svg")
plot(AR_recovery)
dev.off()

#plot(coda::as.mcmc(chain1))

infChain <- data.table::fread(res$history_file,data.table=FALSE)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
#infChain <- data.table::fread(infChainFile,data.table=FALSE)

n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
#beepr::beep(sound=4)
#pdf(paste0(filename,"_infChain.pdf"))
#plot(coda::as.mcmc(n_hist_chain))
#dev.off()

n_hist_chain <- as.data.frame(n_hist_chain[,2:ncol(n_hist_chain)])
hist_chain <- reshape2::melt(n_hist_chain)


alpha <- parTab[parTab$names == "alpha","values"]
beta <- parTab[parTab$names == "beta","values"]

#y=dbb(0:(46*buckets), 46*buckets, alpha, beta)*dbinom(0:(46*buckets),46*buckets,0.5)
#y <- y/sum(y)
y=density_beta_binom(0:(46*buckets), 46*buckets, alpha, beta)
#y <- dbinom(0:(46*buckets),46,0.5)
dat <- data.frame(x=seq(0,46*buckets,by=1),y=y)

#ggplot(hist_chain) + 
#  geom_histogram(aes(x=value,y=..density..),binwidth=1,col="blue",fill="blue",alpha=0.1) + 
#  facet_wrap(~variable) + 
#  theme_bw() +
#  #geom_vline(data=totalInfsDat, aes(xintercept=y),col="red") +
#  geom_line(data=dat,aes(x=x,y=y),col="red",size=0.5)

message(cat("Min ess theta: ", min(effectiveSize(chain1[,which(parTab$fixed==0)+1])),sep="\t"))
message(cat("Min ess infHist: ", min(effectiveSize(n_hist_chain)),sep="\t"))

#message(cat("Min gelman: ", min(effectiveSize(chain1[,which(parTab$fixed==0)+1])),sep="\t"))

infChain <- data.table::fread(res$history_file,data.table=FALSE)
infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
wow <- plyr::ddply(infChain, ~individual, function(x) colSums(x[,!(colnames(x) %in% c("individual","sampno"))])/nrow(x))
wow <- reshape2::melt(wow, id.vars="individual")
wow$variable <- as.integer(wow$variable)
infHist1 <- as.data.frame(infHist)
infHist1$individual <- 1:n_indiv
infHist1 <- infHist1[infHist1$individual %in% 1:10,]

colnames(infHist1) <- c(strainIsolationTimes,"individual")
infHist1 <- reshape2::melt(infHist1, id.vars="individual")
infHist1 <- infHist1[infHist1$value == 1,]
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
realHistProfiles <- as.data.frame(t(apply(infHist, 1, cumsum)))
colnames(realHistProfiles) <- strainIsolationTimes
realHistProfiles$individual <- 1:n_indiv
realHistProfiles <- realHistProfiles[realHistProfiles$individual %in% 1:10,]
realHistProfiles <- reshape2::melt(realHistProfiles,id.vars="individual")
startHistProfiles <- as.data.frame(t(apply(startInf, 1, cumsum)))
colnames(startHistProfiles) <- strainIsolationTimes
startHistProfiles$individual <- 1:n_indiv
startHistProfiles <- startHistProfiles[startHistProfiles$individual %in% 1:10,]
startHistProfiles <- reshape2::melt(startHistProfiles,id.vars="individual")

p1 <- ggplot(quantHist[quantHist$individual %in% seq(1,10,by=1),]) + 
  geom_line(aes(x=as.integer(variable),y=median)) + 
  geom_ribbon(aes(x=as.integer(variable),ymin=lower,ymax=upper), alpha=0.2) + 
  geom_line(data=realHistProfiles[realHistProfiles$individual %in% seq(1,10,by=1),],aes(x=as.integer(variable),y=value), col="blue") +
  geom_line(data=startHistProfiles[startHistProfiles$individual %in% seq(1,10,by=1),],aes(x=as.integer(variable),y=value), col="red") +
  facet_wrap(~individual, scales="free_y")

svg("lambda_sim_200_infectionhistories.svg")
plot(p1)
dev.off()

generate_all_plots(getwd(), mcmcPars["adaptive_period"] + mcmcPars["burnin"], res$chain_file, res$history_file,
                   titreDat, fit_dat, parTab, ages, nIndiv=10,nSamp=100,
                   filename="lambda_sim_200_test")

