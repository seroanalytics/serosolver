library(coda)
library(plyr)

runs <- data.frame("mean"=7,"var"=c(10,25,50,100))
runs <- mefa:::rep.data.frame(runs,2)
runs <- cbind(runs, "buckets"=c(rep(1,4),rep(12,4)))
runs$runName <- paste("var_test",runs$mean,runs$var,runs$buckets,sep="_")
mcmcPars <- c("iterations"=1000000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=1000,"thin"=100,
            "adaptive_period"=500000, "save_block"=500,"thin2"=500,
            "histSampleProb"=1,"switch_sample"=2, "burnin"=0, "moveSize"=5,
            "nInfs"=1,"histProposal"=3,"histOpt"=1)
runs$runName <- as.character(runs$runName)
runs$mean <- as.numeric(runs$mean)
runs$var <- as.numeric(runs$var)
runs$buckets <- as.numeric(runs$buckets)
runs <- runs[,c("runName","mean","var","buckets")]

wds <- unique(as.character(runs$runName))
saveDir <- "~/Documents/Fluscape/serosolver_own/sim_var/"

ess_univ_all <- NULL
max_gelman_univ_all <- NULL
mean_gelman_univ_all <- NULL
mpsrf_gelman_univ_all <- NULL
gelman_hist_all <- NULL
ess_hist_all <- NULL

burnin1 <- 500000
parEstimatesAll <- NULL

index <- 1
for(wd in wds){
  index <- index + 1
  print(wd)
  actualWd <- paste0("~/net/home/serosolver/outputs/",wd)
  buckets <- runs[runs$runName == wd, "buckets"]
  mean <- runs[runs$runName == wd,"mean"]
  var <- runs[runs$runName == wd,"var"]
  
  
  if(length(Sys.glob(file.path(actualWd,"*_chain.pdf"))) < 1){
    max_gelman_univ_all <- c(max_gelman_univ_all, NA)
    mean_gelman_univ_all <- c(mean_gelman_univ_all, NA)
    mpsrf_gelman_univ_all <- c(mpsrf_gelman_univ_all, NA)
    ess_univ_all <- c(ess_univ_all, NA)
  } else {
    chain_files <- Sys.glob(file.path(actualWd, "*_chain.csv"))
    chains_univ <- NULL
    for(i in 1:length(chain_files)){
      tmp <- data.table::fread(chain_files[i],data.table=FALSE)
      vars <- c("sampno","mu", "mu_short", "wane", "tau", "sigma1", "sigma2", "error", 
                "lnlike")
      tmp <- as.mcmc(tmp[tmp$sampno > burnin1,vars])
      chains_univ[[i]] <- tmp
    }
    vars1 <- c("mu", "mu_short", "wane", "tau", "sigma1", "sigma2", "error", 
      "lnlike")
    maxSampno <- min(as.numeric(lapply(chains_univ,function(x) max(x[,"sampno"]))))
    for(i in 1:length(chains_univ)) chains_univ[[i]] <- as.mcmc(chains_univ[[i]][chains_univ[[i]][,"sampno"] < maxSampno,vars1])
      
    chains_univ <- as.mcmc.list(chains_univ)
    
    parEst <- summary(chains_univ)
    lower_quant <- parEst[[2]][,1]
    upper_quant <- parEst[[2]][,5]
    meanPar <- parEst[[1]][,1]
    
    parEstimates <- data.frame("var"=names(meanPar),"mean"=meanPar,"lower"=lower_quant,"upper"=upper_quant,
                               "runName"=wd, "mean1"=mean,"var"=var,"buckets"=buckets)
    parEstimatesAll <- rbind(parEstimatesAll, parEstimates)
    
    chainName1 <- paste0(saveDir, wd,"_univ.pdf")
    
    pdf(chainName1)
    plot(chains_univ)
    dev.off()
    
    ######## Infection histories
    hist_files <- Sys.glob(file.path(actualWd, "*_infectionHistories.csv"))
    hist_chains <- NULL
    for(i in 1:length(hist_files)){
      infChain <- data.table::fread(hist_files[i],data.table=FALSE)
      infChain <- infChain[infChain$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
      infChain <- infChain[infChain$sampno %in% seq(min(infChain$sampno),max(infChain$sampno),by=10),]
      n_infs <- ddply(infChain, ~individual, function(x) summary(rowSums(x[,1:(ncol(x)-2)])))
      n_inf_chain <- ddply(infChain, c("individual","sampno"), function(x) rowSums(x[,1:(ncol(x)-2)]))
      n_hist_chain <- reshape2::dcast(n_inf_chain, sampno~individual, drop=TRUE)
      hist_chains[[i]] <- n_hist_chain
    }
    
    maxSampno <- min(as.numeric(lapply(hist_chains,function(x) max(x[,"sampno"]))))
    for(i in 1:length(hist_chains)) hist_chains[[i]] <- as.mcmc(hist_chains[[i]][hist_chains[[i]][,"sampno"] < maxSampno,2:ncol(hist_chains[[i]])])
    hist_chains <- as.mcmc.list(hist_chains)
    
    chainName2 <- paste0(saveDir, wd,"_hist.pdf")
    
    pdf(chainName2)
    plot(hist_chains)
    dev.off()
    ess_hist <- min(effectiveSize(hist_chains))
    gelman_hist <- gelman.diag(hist_chains)
    
    gelman_hist_all <- c(gelman_hist_all, max(gelman_hist$psrf[,2]))
    ess_hist_all <- c(ess_hist_all, ess_hist)
    
    ess_univ <- min(effectiveSize(chains_univ))
    gelman_univ <- gelman.diag(chains_univ)
    ess_univ_all <- c(ess_univ_all, ess_univ)
    max_gelman_univ <- max(gelman_univ$psrf[,2])
    mean_gelman_univ <- mean(gelman_univ$psrf[,2])
    mpsrf_gelman_univ <- gelman_univ$mpsrf
    max_gelman_univ_all <- c(max_gelman_univ_all, max_gelman_univ)
    mean_gelman_univ_all <- c(mean_gelman_univ_all, mean_gelman_univ)
    mpsrf_gelman_univ_all <- c(mpsrf_gelman_univ_all, mpsrf_gelman_univ)
  }
}
parEstimatesAll$version <- as.factor(parEstimatesAll$version)
parEstimatesAll$histProposal <- as.factor(parEstimatesAll$histProposal)
parEstimatesAll$histOpt <- as.factor(parEstimatesAll$histOpt)
parTab <- read.csv("~/net/home/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab <- parTab[parTab$fixed == 0,]
colnames(parTab)[1] <- "var"
p <- ggplot(parEstimatesAll) + geom_pointrange(aes(x=runName, y = mean,ymin=lower,ymax=upper, col=version,shape=histProposal)) + 
  geom_hline(data=parTab, aes(yintercept=values),linetype="dashed")+
  facet_wrap(~var,scales="free_y") + 
  theme(axis.text.x=element_text(angle=90,hjust=0.5,size=6))

pdf("~/Documents/Fluscape/serosolver_own/sim_var/simrecover.pdf")
plot(p)
dev.off()
svg("~/Documents/Fluscape/serosolver_own/sim_var/simrecover.svg", width=10,height=9)
plot(p)
dev.off()

res <- data.frame(wds, runs$runName, runs$version, runs$histProposal, runs$histOpt,ess_univ_all,           
                  max_gelman_univ_all, mean_gelman_univ_all, mpsrf_gelman_univ_all,
                  gelman_hist_all, ess_hist_all)

write.table(res,"~/Documents/Fluscape/serosolver_own/sim_priors/simrecover_res.csv",row.names=FALSE,sep=",")
write.table(parEstimatesAll,"~/Documents/Fluscape/serosolver_own/sim_priors/simrecover_parEst.csv",row.names=FALSE,sep=",")
