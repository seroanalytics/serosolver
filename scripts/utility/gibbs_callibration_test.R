library(coda)
library(data.table)
library(ggplot2)
library(reshape2)
library(plyr)
setwd("/media/james/F40474B404747B84/Storage/serosolver/outputs_15082018/to_analyse1")
topdir <- getwd()
dirs <- list.dirs(path=topdir,recursive=FALSE,full.names=FALSE)
tmps <- NULL
for(dir in dirs){
  setwd(topdir)
  setwd(dir)
  files <- list.files(pattern="_chain.csv")
  tmp <- NULL
  for(file in files){
    chain <- data.table::fread(file,data.table=FALSE)
    chain <- chain[chain$sampno < 240000 & chain$sampno > 50000 ,c("mu","mu_short","tau","sigma1","sigma2","error","wane")]
    tmp[[file]] <- as.mcmc(chain)
  }
  pdf(paste0(topdir,"/",dir,"_chains.pdf"))
  plot(as.mcmc.list(tmp))
  dev.off()
  tmps[[dir]] <- tmp
}
#lapply(tmps, function(x) lapply(x, function(y) range(y[,"sampno"])))
res <- lapply(tmps, function(x) summary(effectiveSize(as.mcmc.list(x))))
res <- do.call("rbind",res)
res <- cbind("sampPropn"=c(0.5,0.5,0.5,0.5,1,1,1,1),"histPropn"=c(0.5,0.5,1,1,0.5,0.5,1,1),"moveSize"=c(2,5,2,5,2,5,2,5),res)
res_chain <- res
tmps <- NULL
for(dir in dirs){
  setwd(topdir)
  setwd(dir)
  files <- list.files(pattern="infectionHistories.csv")
  tmp <- NULL
  print(dir)
  index <- 1
  print(files)
  for(file in files){
    print(index)
    infChain <- data.table::fread(file)
    infChain <- infChain[infChain$sampno < 240000 & infChain$sampno > 50000 ,]
    n_strain <- max(infChain$j)
    data.table::setkey(infChain, "j","sampno")
    n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
    n_inf_chain <- cbind(n_inf_chain,chain=index)
    index <- index + 1
    tmp <- rbind(tmp, n_inf_chain)
    #tmp[[file]] <- as.mcmc(chain)
  }
  tmp$chain <- as.factor(tmp$chain)
  pdf(paste0(topdir,"/",dir,"_inf_chains.pdf"))
  p1 <- ggplot(tmp) + geom_line(aes(x=sampno,y=V1,col=chain)) + facet_wrap(~j)
  print(p1)
  dev.off()
  ess <- t(ddply(tmp,~j, function(x) effectiveSize(x$V1)))
  tmps <- rbind(tmps,ess[2,])
}
res <- tmps
res <- cbind("sampPropn"=c(0.5,0.5,0.5,0.5,1,1,1,1),"histPropn"=c(0.5,0.5,1,1,0.5,0.5,1,1),"moveSize"=c(2,5,2,5,2,5,2,5),res)

res <- data.frame(res)
res$run <- dirs
res2 <- melt(res,id.vars=c("sampPropn","histPropn","moveSize","run"))
res2$sampPropn <- as.factor(res2$sampPropn)
res2$histPropn <- as.factor(res2$histPropn)
res2$moveSize <- as.factor(res2$moveSize)
res2$run <- as.factor(res2$run)

p1 <- ggplot(res2) + geom_bar(aes(x=variable,y=value,fill=moveSize),stat="identity",position="dodge") + facet_wrap(~run)
p2 <- ggplot(res2) + geom_bar(aes(x=variable,y=value,fill=histPropn),stat="identity",position="dodge") + facet_wrap(~run)
p3 <- ggplot(res2) + geom_bar(aes(x=variable,y=value,fill=sampPropn),stat="identity",position="dodge") + facet_wrap(~run)

cowplot::plot_grid(p1,p2,p3,ncol=3)

