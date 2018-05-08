library(reshape2)
library(plyr)
library(ggplot2)
source("model_funcs.R")
source("probability_funcs.R")
source("proposal_funcs.R")
source("mcmc_funcs.R")

## Input parameters
n <- 20
coin_probs <- runif(n)
indivs <- 10
samps <- seq(1,n, by=2)
pars <- c(4, 0.3, 1, coin_probs)
fixed <- c(0,0,1,rep(0,n))
covMat <- diag(length(fixed[which(fixed==0)]))*0.01
iter <- 10000


## Setup parameter names and simulated data
parNames <- c("boost","sigma","error",paste0("coin.",1:n))
coin_results <- sapply(coin_probs, function(x) sample(c(0,1),indivs,prob=c(1-x,x),replace=TRUE))
data_suggested_coins <- colSums(coin_results)/nrow(coin_results)
dat <- coin_toss_group(pars, coin_results)
dat <- measurement_error_group(pars,dat)


print(paste0("Coin probs: ", paste0(signif(coin_probs,3),collapse=" ")))
print(paste0("Data suggest coin values: ", paste0(data_suggested_coins,collapse=" ")))

## Starting conditions and run MCMC
startPars <- pars
start_coins <- matrix(sample(c(0,1),n*indivs,replace=TRUE),nrow=indivs)
res <- run_MCMC_group(pars, fixed, dat, start_coins,iter, covMat,thin=10, samps)

## Look at MCMC chain of process parameters
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
  geom_point(data=real_coins, aes(x=variable,y=y),col="red") +
  geom_point(data=real_coins_data,aes(x=variable,y=y),col="purple") +
  scale_y_continuous(limits=c(0,1))



liks <- NULL
liks1 <- NULL
pars1 <- pars
priors <- NULL
for(i in 1:100){
  pars1[4] <- (i-1)/100
  priors[i] <- prior(pars1)
  liks[i] <- posterior_group(pars1, coin_results, dat)
  liks1[i] <- dbinom(sum(coin_results[,1]), nrow(coin_results), pars1[4],log=TRUE)
}
#plot(liks1)
#lines(liks+log(choose(indivs, sum(coin_results[,2]))),col="red")

