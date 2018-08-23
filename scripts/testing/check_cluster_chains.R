library(coda)
chains <- list.files(pattern="chain.csv")
omg <- NULL
for(chain in chains){
  tmp <- read.csv(chain)
  tmp <- tmp[tmp$sampno > 10000 & tmp$sampno < 19900,]
  #c("mu","mu_short","tau","wane","sigma1","sigma2","error","lnlike")
  #tmp <- tmp[tmp$sampno < 8000,c("mu","mu_short","tau","wane","sigma1","sigma2","error","lnlike")]
  omg[[chain]] <- as.mcmc(tmp)
}
omg <- as.mcmc.list(omg)
#omg <- as.mcmc.list(omg[c(1,2,5,6,8)])
#omg <- omg[-1]
pdf("tmp.pdf")
plot(omg)
dev.off()

library(ggplot2)
library(data.table)
library(plyr)
n_indiv <- 1000
indivs <- sample(n_indiv, 10)
sampd <- sample(n_indiv,20)
all_chain <- NULL
infHistFiles <- list.files(pattern="infectionHistories.csv")
for(i in 1:length(infHistFiles)){
  infChainFile <- infHistFiles[i]
  infChain <- data.table::fread(infChainFile)
   infChain <- infChain[infChain$sampno > 50000 & infChain$sampno < 250000,]
   n_strain <- max(infChain$j)
   data.table::setkey(infChain, "j","sampno")
   n_inf_chain <- infChain[,list(V1=sum(x)),by=key(infChain)]
   n_inf_chain <- cbind(n_inf_chain,"chain"=i)
   all_chain <- rbind(all_chain, n_inf_chain)
} 
#all_chain <- all_chain[all_chain$chain %in% c(1,2,5,6,8),]
all_chain$chain <- as.factor(all_chain$chain)
#all_chain <- all_chain[all_chain$chain != 1,]
inf_chain_p <- ggplot(all_chain) + geom_line(aes(x=sampno,y=V1,col=chain)) + facet_wrap(~j)
devtools::load_all("~/Documents/Fluscape/serosolver")
titreFiles <- list.files(pattern="titreDat")
titreDat <- read.csv(titreFiles[1])
ageFiles <- list.files(pattern="ages")
ages <- read.csv(ageFiles[1])
AR_p <- plot_attack_rates(infChain, titreDat, ages, seq(1968, 2015, by=1)) 

library(coda)
chains <- list.files(pattern="multivariate_chain.csv")
omg <- NULL
for(chain in chains){
  tmp <- read.csv(chain)
  tmp <- tmp[tmp$sampno > 50000,]
  omg <- rbind(omg,tmp)
  
}
