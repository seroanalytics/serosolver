library(coda)
runs <- read.csv("~/net/home/serosolver/inputs/runs_18012018.csv",stringsAsFactors=FALSE)

runs <- runs[runs$buckets %in% c(1,4,12),]
wds <- as.character(runs$runName)
saveDir <- "~/Documents/Fluscape/serosolver_own/sim1/"

ess_univ_all <- NULL
ess_multi_all <- NULL
max_gelman_univ_all <- NULL
mean_gelman_univ_all <- NULL
mpsrf_gelman_univ_all <- NULL
max_gelman_multi_all <- NULL
mean_gelman_multi_all <- NULL
mpsrf_gelman_multi_all <- NULL

burnin1 <- burnin2 <- 300000

parEstimatesAll <- NULL

index <- 1
for(wd in wds){
  index <- index + 1
  print(wd)
  actualWd <- paste0("~/net/home/serosolver/outputs/",wd)
  buckets <- runs[runs$runName == wd,"buckets"]
  alpha <- runs[runs$runName == wd,"alpha"]
  histSampProb <- runs[runs$runName == wd,"histSampleProb"]
  
  if(length(Sys.glob(file.path(actualWd,"*_theta_chain.pdf"))) != 3){
    max_gelman_univ_all <- c(max_gelman_univ_all, NA)
    mean_gelman_univ_all <- c(mean_gelman_univ_all, NA)
    mpsrf_gelman_univ_all <- c(mpsrf_gelman_univ_all, NA)
    
    max_gelman_multi_all <- c(max_gelman_multi_all, NA)
    mean_gelman_multi_all <- c(mean_gelman_multi_all, NA)
    mpsrf_gelman_multi_all <- c(mpsrf_gelman_multi_all, NA)
    ess_univ_all <- c(ess_univ_all, NA)
    ess_multi_all <- c(ess_multi_all, NA)
  } else {
    chain_files_univ <- Sys.glob(file.path(actualWd, "*_univariate_chain.csv"))
    chain_files_multi <- Sys.glob(file.path(actualWd, "*_multivariate_chain.csv"))
  
    chains_univ <- NULL
    for(i in 1:length(chain_files_univ)){
      tmp <- data.table::fread(chain_files_univ[i],data.table=FALSE)
      vars <- c("sampno","mu", "mu_short", "wane", "tau", "sigma1", "sigma2", "error", 
                "lnlike")
      tmp <- as.mcmc(tmp[tmp$sampno > burnin1,vars])
      chains_univ[[i]] <- tmp
    }
    vars1 <- c("mu", "mu_short", "wane", "tau", "sigma1", "sigma2", "error", 
      "lnlike")
    maxSampno <- min(as.numeric(lapply(chains_univ,function(x) max(x[,"sampno"]))))
    for(i in 1:length(chains_univ)) chains_univ[[i]] <- as.mcmc(chains_univ[[i]][chains_univ[[i]][,"sampno"] < maxSampno,vars1])
      
    
    chains_multi <- NULL
    for(i in 1:length(chain_files_multi)){
      tmp <- data.table::fread(chain_files_multi[i],data.table=FALSE)
      vars <- c("sampno","mu", "mu_short", "wane", "tau", "sigma1", "sigma2", "error", 
                "lnlike")
      tmp <- as.mcmc(tmp[tmp$sampno > burnin2,vars])
      chains_multi[[i]] <- tmp
    }
    maxSampno <- min(as.numeric(lapply(chains_multi,function(x) max(x[,"sampno"]))))
    for(i in 1:length(chains_multi)) chains_multi[[i]] <- as.mcmc(chains_multi[[i]][chains_multi[[i]][,"sampno"] < maxSampno,vars1])
    
    
    chains_multi <- as.mcmc.list(chains_multi)
    chains_univ <- as.mcmc.list(chains_univ)
    
    parEst <- summary(chains_multi)
    lower_quant <- parEst[[2]][,1]
    upper_quant <- parEst[[2]][,5]
    meanPar <- parEst[[2]][,1]
    
    parEstimates <- data.frame("var"=names(meanPar),"mean"=meanPar,"lower"=lower_quant,"upper"=upper_quant,
                               "runName"=wd,"buckets"=buckets,"alpha"=alpha,"histSampProb"=histSampProb)
    parEstimatesAll <- rbind(parEstimatesAll, parEstimates)
    
    chainName1 <- paste0(saveDir, wd,"_univ.pdf")
    chainName2 <- paste0(saveDir, wd,"_multi.pdf")
    
    pdf(chainName1)
    plot(chains_univ)
    dev.off()
    
    pdf(chainName2)
    plot(chains_multi)
    dev.off()
    
    ess_univ <- min(effectiveSize(chains_univ))
    gelman_univ <- gelman.diag(chains_univ)
    
    ess_multi <- min(effectiveSize(chains_multi))
    gelman_multi <- gelman.diag(chains_multi)
    
    ess_univ_all <- c(ess_univ_all, ess_univ)
    ess_multi_all <- c(ess_multi_all, ess_multi)
    
    max_gelman_univ <- max(gelman_univ$psrf[,2])
    mean_gelman_univ <- mean(gelman_univ$psrf[,2])
    mpsrf_gelman_univ <- gelman_univ$mpsrf
    
    max_gelman_multi <- max(gelman_multi$psrf[,2])
    mean_gelman_multi <- mean(gelman_multi$psrf[,2])
    mpsrf_gelman_multi <- gelman_multi$mpsrf
    max_gelman_univ_all <- c(max_gelman_univ_all, max_gelman_univ)
    mean_gelman_univ_all <- c(mean_gelman_univ_all, mean_gelman_univ)
    mpsrf_gelman_univ_all <- c(mpsrf_gelman_univ_all, mpsrf_gelman_univ)
    
    max_gelman_multi_all <- c(max_gelman_multi_all, max_gelman_multi)
    mean_gelman_multi_all <- c(mean_gelman_multi_all, mean_gelman_multi)
    mpsrf_gelman_multi_all <- c(mpsrf_gelman_multi_all, mpsrf_gelman_multi)
  }
  
}
parEstimatesAll$alpha <- as.factor(parEstimatesAll$alpha)
parEstimatesAll$buckets <- as.factor(parEstimatesAll$buckets)
parTab <- read.csv("~/net/home/serosolver/inputs/parTab.csv",stringsAsFactors=FALSE)
parTab <- parTab[parTab$fixed == 0,]
colnames(parTab)[1] <- "var"
p <- ggplot(parEstimatesAll) + geom_pointrange(aes(x=runName, y = mean,ymin=lower,ymax=upper, col=buckets,shape=alpha)) + 
  geom_hline(data=parTab, aes(yintercept=values),linetype="dashed")+
  facet_wrap(~var,scales="free_y") + 
  theme(axis.text.x=element_text(angle=90,hjust=0.5,size=6))

pdf("~/Documents/Fluscape/serosolver_own/simrecover.pdf")
plot(p)
dev.off()
svg("~/Documents/Fluscape/serosolver_own/simrecover.svg", width=10,height=9)
plot(p)
dev.off()


res <- data.frame(wds, runs$N, runs$histSampleProb, runs$switch_sample, runs$buckets,runs$alpha,runs$beta,ess_univ_all, 
                  ess_multi_all, 
             max_gelman_univ_all, mean_gelman_univ_all, mpsrf_gelman_univ_all,
             max_gelman_multi_all, mean_gelman_multi_all, mpsrf_gelman_multi_all)
write.table(res,"~/Documents/Fluscape/serosolver_own/simrecover_res.csv",row.names=FALSE,sep=",")
write.table(parEstimatesAll,"~/Documents/Fluscape/serosolver_own/simrecover_parEst.csv",row.names=FALSE,sep=",")
