devtools::load_all("~/Documents/Fluscape/serosolver")
source("~/Documents/Fluscape/flu-model/sero_model/sero_functions.R")
library(reshape2)
library(plyr)
library(ggplot2)


#tmp <- inf_hist_prop_cpp(as.matrix(infHist), sampled, ageMask, moveSizes, nInfs, 100, 100,rep(0.2,10))



n_indiv <- 25
infHist <- infHist1 <- infHist2 <- infHist3 <- infHist4 <- matrix(0, nrow=n_indiv, ncol=10,byrow=TRUE)
sampled <- 1:n_indiv
strainIsolationTimes <- 1:10
ageMask <- sample(1:5,n_indiv,replace=TRUE)
ageMask <- rep(1, n_indiv)
nInfs <- rep(1,n_indiv)
moveSizes <- rep(1, n_indiv)

omgs <- omgs1 <- omgs2 <- omgs3 <- omgs4 <- list()

n <- 10000
infHists <- infHists1 <- infHists2 <- infHists3 <- infHists4 <- matrix(ncol=13,nrow=n_indiv*n)

for(i in 1:n){
  if(i%%100 == 0) print(i)
  infHist1 <- as.data.frame(infection_history_betabinom_group(infHist1, sampled, ageMask, moveSizes, nInfs, 1, 1))
  infHists1[(((i-1)*n_indiv)+1):(i*n_indiv),1:10] <- as.matrix(infHist1)
  infHists1[(((i-1)*n_indiv)+1):(i*n_indiv),11] <- 1:n_indiv
  infHists1[(((i-1)*n_indiv)+1):(i*n_indiv),12] <- i 
  tmp1 <- data.frame(infs=rowSums(infHist1),indiv=1:n_indiv,sample=i)
  omgs1[[i]] <- tmp1
  
  infHist2 <- as.data.frame(infection_history_proposal(infHist2, sampled, strainIsolationTimes, ageMask, nInfs)[[1]])
  infHists2[(((i-1)*n_indiv)+1):(i*n_indiv),1:10] <-  as.matrix(infHist2)
  infHists2[(((i-1)*n_indiv)+1):(i*n_indiv),11] <- 1:n_indiv
  infHists2[(((i-1)*n_indiv)+1):(i*n_indiv),12] <- i 
  tmp2 <- data.frame(infs=rowSums(infHist2),indiv=1:n_indiv,sample=i)
  omgs2[[i]] <- tmp2
  
  infHist <- as.data.frame(SampleHistory(infHist, sampled, 10, ageMask, strainIsolationTimes,rep(1,n_indiv)))
  infHists[(((i-1)*n_indiv)+1):(i*n_indiv),1:10] <-  as.matrix(infHist)
  infHists[(((i-1)*n_indiv)+1):(i*n_indiv),11] <- 1:n_indiv
  infHists[(((i-1)*n_indiv)+1):(i*n_indiv),12] <- i 
  tmp <- data.frame(infs=rowSums(infHist),indiv=1:n_indiv,sample=i)
  omgs[[i]] <- tmp
  
  infHist3 <- as.data.frame(inf_hist_prop_cpp(as.matrix(infHist3), sampled, ageMask, moveSizes, nInfs, 1, 1,runif(n_indiv)))
  infHists3[(((i-1)*n_indiv)+1):(i*n_indiv),1:10] <-  as.matrix(infHist3)
  infHists3[(((i-1)*n_indiv)+1):(i*n_indiv),11] <- 1:n_indiv
  infHists3[(((i-1)*n_indiv)+1):(i*n_indiv),12] <- i 
  tmp3 <- data.frame(infs=rowSums(infHist3),indiv=1:n_indiv,sample=i)
  omgs3[[i]] <- tmp3
  
  infHist4 <- as.data.frame(infection_history_symmetric(as.matrix(infHist4), sampled, ageMask, moveSizes, nInfs, runif(n_indiv)))
  infHists4[(((i-1)*n_indiv)+1):(i*n_indiv),1:10] <-  as.matrix(infHist4)
  infHists4[(((i-1)*n_indiv)+1):(i*n_indiv),11] <- 1:n_indiv
  infHists4[(((i-1)*n_indiv)+1):(i*n_indiv),12] <- i 
  tmp4 <- data.frame(infs=rowSums(infHist4),indiv=1:n_indiv,sample=i)
  omgs4[[i]] <- tmp4
}
omgs <- do.call("rbind",omgs)
omgs1 <- do.call("rbind",omgs1)
omgs2 <- do.call("rbind",omgs2)
omgs3 <- do.call("rbind",omgs3)
omgs4 <- do.call("rbind",omgs4)


infHists <- as.data.frame(infHists)
infHists1 <- as.data.frame(infHists1)
infHists2 <- as.data.frame(infHists2)
infHists3 <- as.data.frame(infHists3)
infHists4 <- as.data.frame(infHists4)
infHists[,13] <- "AK"
infHists1[,13] <- "beta_binom"
infHists2[,13] <- "AK_mine"
infHists3[,13] <- "BB_cpp"
infHists4[,13] <- "symmetric"
colnames(infHists) <- colnames(infHists1) <-colnames(infHists2) <-colnames(infHists3) <-colnames(infHists4) <- c(1:10, "indiv","samp","ver")

infHist5 <- matrix(0, nrow=n_indiv, ncol=10,byrow=TRUE)
infHists5 <- matrix(ncol=13,nrow=n_indiv*n)
omgs5 <- NULL
#ageMask <- rep(5,n_indiv)
#infHist5[,1:5] <- 0
for(i in 1:n){
  if(i%%100 == 0) print(i)
  infHist5 <- as.data.frame(inf_hist_prop_cpp(as.matrix(infHist5), sampled, ageMask, moveSizes, nInfs, 1, 1,runif(n_indiv, 0, 0.5)))
  #infHist5 <- as.data.frame(infection_history_betabinom_group(infHist5, sampled, ageMask, moveSizes, nInfs, 1, 1))
  #infHist5 <- as.data.frame(infection_history_symmetric(as.matrix(infHist5), sampled, ageMask, moveSizes, nInfs, runif(n_indiv)))
  infHists5[(((i-1)*n_indiv)+1):(i*n_indiv),1:10] <-  as.matrix(infHist5)
  infHists5[(((i-1)*n_indiv)+1):(i*n_indiv),11] <- 1:n_indiv
  infHists5[(((i-1)*n_indiv)+1):(i*n_indiv),12] <- i 
  tmp5 <- data.frame(infs=rowSums(infHist5),indiv=1:n_indiv,sample=i)
  omgs5 <- rbind(omgs5, tmp5)
}
infHists5 <- as.data.frame(infHists5)
colnames(infHists5) <- c(1:10, "indiv","samp","ver")
infHists5[,13] <- "binom_bb"
omgs5$ver <- "binom_bb"

infHists_comb <- as.data.frame(rbind(infHists,infHists1, infHists2, infHists3, infHists4, infHists5))
#infHists_comb <- infHists5[complete.cases(infHists5),]

greb <- melt(infHists_comb, id.vars=c("ver","samp","indiv"))
AR <- ddply(greb, .(variable, ver, samp), function(x) sum(x$value)/n_indiv)
AR1 <- ddply(greb, .(variable, ver, samp), function(x) sum(x$value))
res <- ddply(AR1, .(variable,ver), function(x) quantile(x$V1,c(0.025,0.5,0.975)))
res <- res[order(res$ver,res$variable),]
n_alive <- sapply(strainIsolationTimes, function(x) length(ageMask[ageMask <=x]))
#n_alive <- 25
res[,3] <- res[,3]/rep(n_alive,1)
res[,4] <- res[,4]/rep(n_alive,1)
res[,5] <- res[,5]/rep(n_alive,1)
ggplot(res) + geom_pointrange(aes(x=variable,y=`50%`,ymax=`97.5%`,ymin=`2.5%`)) + facet_wrap(~ver) + scale_y_continuous(limits=c(0,1))

for(var in unique(AR1$variable)){
  AR1[AR1$variable == var, "V1"] <- AR1[AR1$variable == var, "V1"]/n_alive[as.numeric(var)]
}
AR_densities <- ggplot(AR1) + 
  geom_density(aes(x=V1)) + 
  #scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2),labels=seq(0,1,by=0.2)) +
  facet_grid(ver~variable) + 
  theme_bw()

omgs$ver <- "AK"
omgs1$ver <- "beta_binom"
omgs2$ver <- "AK_mine"
omgs3$ver <- "BB_cpp"
omgs4$ver <- "symmetric"
wow <- rbind(omgs, omgs1, omgs2, omgs3, omgs4, omgs5)
wow$infs <- as.integer(wow$infs)
#wow <- omgs5
#wow[wow$indiv %in% which(ageMask >= 5),]
indiv_total <- ggplot(wow) + geom_histogram(aes(x=infs,fill=ver),binwidth=1) + 
    facet_grid(ver~indiv) + 
  #facet_wrap(~indiv)+
  scale_x_continuous(expand=c(0,0),
                     limits=c(0,10),labels=seq(0,10,by=2),breaks=seq(0,10,by=2))










