library(coda)
setwd("/media/james/Storage 2/serosolver_21Aug2018/callibrate")
dirs <- list.dirs(recursive = FALSE,full.names = FALSE)
#dirs <- list.files(pattern="sim")
library(data.table)
topwd <- getwd()
ess_all <- NULL
gelman_all <- NULL
ess_hist <- NULL
times_all <- NULL
for(dir in dirs){
  setwd(topwd)
  setwd(dir)
  files <- list.files(pattern="chain.csv")  
  infHistFiles <- list.files(pattern="infectionHistories.csv")
  chains <- NULL
  all_chain <- NULL
  times <- NULL
  for(i in 1:length(files)){
    tmp <- fread(files[i],data.table=FALSE)
    tmp <- as.mcmc(tmp[tmp$sampno > 50000,c("sampno","mu","mu_short","tau","sigma1","sigma2","error")])
    chains[[i]] <- tmp
    infChainFile <- infHistFiles[i]
    infChain <- data.table::fread(infChainFile)
    infChain <- infChain[infChain$sampno > 50000,]
    n_strain <- max(infChain$j)
    data.table::setkey(infChain, "j","sampno")
    n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
    n_inf_chain <- cbind(n_inf_chain,"chain"=i)
    all_chain <- rbind(all_chain, n_inf_chain)
    times[[i]] <- as.numeric(strsplit(as.character(read.delim(list.files(pattern="time")[1])[1,1]),split=" ")[[1]][4])
  }
  times_all[[dir]] <- mean(times)
  print(paste0("dir: ", dir, " is ", mean(times)))
  ess_hist[[dir]] <- t(plyr::ddply(all_chain, ~j, function(x) effectiveSize(x$V1)))[2,]
  chains <- as.mcmc.list(chains)
  gelman_all[[dir]] <- gelman.diag(chains)
  ess_all[[dir]] <- effectiveSize(chains)
  setwd(topwd)
  pdf(paste0(dir,"_chains.pdf"))
  plot(chains)
  dev.off()
}

ess_sigma_comb <- as.data.frame(do.call("rbind",ess_all))
ess_hist_comb <- as.data.frame(do.call("rbind",ess_hist))
ess_hist_comb$runName <- row.names(ess_hist_comb)
times_all <- c(times_all)

all_res <- cbind(ess_hist_comb,ess_sigma_comb,times_all)
all_res[all_res$times_all < 5,"times_all"] <- all_res[all_res$times_all < 5,"times_all"]*24

tmp <- reshape2::melt(ess_hist_comb)
library(ggplot2)
tmp$variable <- as.integer(tmp$variable)

ggplot(tmp) + geom_line(aes(x=variable,y=value,col=runName))
