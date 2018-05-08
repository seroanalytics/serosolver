source("model_funcs.R")
source("probability_funcs.R")
source("proposal_funcs.R")
source("mcmc_funcs.R")

library(reshape2)
library(plyr)
library(ggplot2)

n <- 10
indivs <- 20
coin_probs1 <- runif(n)
coin_probs <- matrix(sapply(coin_probs1, function(x) rnorm(indivs,x,0.1)),nrow=indivs)
coin_probs[coin_probs < 0] <- 0
coin_probs[coin_probs > 1] <- 1
pars <- c(4, 0.3, 1)
coin_results <- t(apply(coin_probs, 1, function(x) sapply(x, function(y) sample(c(0,1),1,prob=c(1-y,y),replace=TRUE))))

parNames <- c("boost","sigma","error")
data_suggested_coins <- colSums(coin_results)/nrow(coin_results)
#print(paste0("Coin probs: ", paste0(signif(coin_probs,3),collapse=" ")))
print(paste0("Data suggest coin values: ", paste0(data_suggested_coins,collapse=" ")))
dat <- coin_toss_group(pars, coin_results)
dat <- measurement_error_group(pars,dat)
fixed <- c(0,0,1)
covMat <- diag(length(fixed[which(fixed==0)]))*0.01
iter <- 10000
start_coins <- matrix(sample(c(0,1),n*indivs,replace=TRUE),nrow=indivs)
res <- run_MCMC_group(pars, fixed, dat, start_coins,coin_probs, iter, covMat,thin=10)
chain <- res[[2]]
chain <- chain[,c(1,which(fixed == 0)+1)]
colnames(chain) <- c("sampno",parNames[which(fixed==0)])
chain <- chain[chain[,"sampno"] > 0.1*iter,]
plot(coda::as.mcmc(chain))
chain <- melt(chain)
colnames(chain) <- c("x","parameter","value")
real_pars <- data.frame(parameter=parNames[which(fixed==0)],y=pars[which(fixed==0)])
p1 <- ggplot(chain) + 
  geom_density(aes(x=value,y=..scaled..)) +
  geom_vline(data=real_pars,aes(xintercept=y),col="red",linetype="dashed")+
  facet_wrap(~parameter, scales="free_x")

coin_chain <- as.data.frame(res[[3]])
colnames(coin_chain) <- c("sampno",1:n,"indiv")
coin_chain <- coin_chain[coin_chain$sampno > 0.1*iter,]
allRes <- NULL
for(j in 1:n){
  print(j)
  x <- NULL
  i <- 1
  for(sampno in unique(coin_chain$sampno)){
    tmpChain <- coin_chain[coin_chain$sampno == sampno,]
    x[i] <- sum(tmpChain[,colnames(tmpChain) == j])/indivs
    i <- i+ 1
  }
  final <- quantile(x,c(0.025,0.5,0.975))
  allRes <- rbind(allRes, final)
}

tmp <- melt(coin_chain[,1:(ncol(coin_chain)-1)], id.vars="sampno")
tmp <- ddply(tmp, .(sampno, variable), function(x) sum(x$value)/indivs)
tmp1 <- dcast(tmp, sampno~variable)
tmp <- ddply(tmp, ~variable, function(x) quantile(x$V1, c(0.025,0.5,0.975)))

inferred_coins <- plyr::ddply(coin_chain, ~indiv, colMeans)
inferred_coins <- inferred_coins[,2:ncol(inferred_coins)]
wow <- melt(inferred_coins,id.vars="indiv")
omg <- melt(coin_results)
colnames(omg) <- c("indiv","variable","value")
omg <- omg[omg$value == 1,]
p2 <- ggplot(wow[wow$indiv %in% 1:10,]) + 
  geom_point(aes(x=variable,y=value)) + 
  geom_vline(data=omg[omg$indiv %in% 1:10,],aes(xintercept=variable),col="red",linetype="dashed",alpha=0.4) +
  facet_wrap(~indiv)

real_coins <- data.frame(variable=1:n, y=coin_probs)
real_coins_data <- data.frame(variable=1:n, y=colSums(coin_results)/nrow(coin_results))
tmp2 <- chain[chain$parameter %in% paste0("coin.",1:n),]
colnames(tmp2) <- c("x","variable","value")
tmp2$variable <- factor(tmp2$variable, levels=unique(as.character(tmp2$variable)))
tmp2$variable <- as.integer(tmp2$variable)
tmp2 <- ddply(tmp2, ~variable, function(x) quantile(x$value, c(0.025,0.5,0.975)))
tmp$variable <- as.integer(tmp$variable)
p3 <- ggplot(tmp) + 
  geom_pointrange(aes(x=variable-0.2,ymin=`2.5%`,y=`50%`,ymax=`97.5%`)) + 
  #geom_pointrange(data=tmp2,aes(x=variable+0.2,ymin=`2.5%`,y=`50%`,ymax=`97.5%`),col="blue") + 
  #geom_point(data=real_coins, aes(x=variable,y=y),col="red") +
  geom_point(data=real_coins_data,aes(x=variable,y=y),col="purple") +
  scale_y_continuous(limits=c(0,1))







coin_prob_chain <- as.data.frame(res[[4]])
colnames(coin_prob_chain) <- c("sampno",1:n,"indiv")
#coin_prob_chain <- coin_prob_chain[coin_prob_chain$sampno > 0.1*iter,]
coin_prob_chain <- melt(coin_prob_chain,id.vars=c("sampno","indiv"))
ggplot(coin_prob_chain[coin_prob_chain$indiv==1,]) + geom_line(aes(x=sampno,y=value)) + facet_wrap(~variable)
tmp1 <- dcast(coin_prob_chain, sampno~variable,sum)
tmp <- melt(coin_prob_chain[,1:(ncol(coin_prob_chain)-1)], id.vars="sampno")
#tmp <- ddply(tmp, .(sampno, variable), function(x) sum(x$value)/indivs)
tmp1 <- dcast(tmp, sampno~variable, sum)
tmp <- ddply(tmp, ~variable, function(x) quantile(x$V1, c(0.025,0.5,0.975)))

inferred_coins <- plyr::ddply(coin_chain, ~indiv, colMeans)
inferred_coins <- inferred_coins[,2:ncol(inferred_coins)]
wow <- melt(inferred_coins,id.vars="indiv")
omg <- melt(coin_results)
colnames(omg) <- c("indiv","variable","value")
omg <- omg[omg$value == 1,]
p2 <- ggplot(wow[wow$indiv %in% 1:10,]) + 
  geom_point(aes(x=variable,y=value)) + 
  geom_vline(data=omg[omg$indiv %in% 1:10,],aes(xintercept=variable),col="red",linetype="dashed",alpha=0.4) +
  facet_wrap(~indiv)

