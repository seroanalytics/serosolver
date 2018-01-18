library(coda)
runs <- read.csv("~/net/home/serosolver/inputs/runs_18012018.csv",stringsAsFactors=FALSE)
wds <- as.character(runs$runName)

saveDir <- "~/Documents/Fluscape/serosolver_own/sim/"

ess_univ_all <- NULL
ess_multi_all <- NULL
max_gelman_univ_all <- NULL
mean_gelman_univ_all <- NULL
mpsrf_gelman_univ_all <- NULL
max_gelman_multi_all <- NULL
mean_gelman_multi_all <- NULL
mpsrf_gelman_multi_all <- NULL

burnin1 <- burning2 <- 550000
for(wd in wds){
  print(wd)
  actualWd <- paste0("~/net/home/serosolver/outputs/",wd)
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

res <- data.frame(wds, runs$N, runs$histSampleProb, runs$switch_sample, runs$buckets,runs$alpha,runs$beta,ess_univ_all, 
                  ess_multi_all, 
             max_gelman_univ_all, mean_gelman_univ_all, mpsrf_gelman_univ_all,
             max_gelman_multi_all, mean_gelman_multi_all, mpsrf_gelman_multi_all
             )
