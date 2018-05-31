library(reshape2)
library(plyr)
library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(5)  #change the 2 to your number of CPU cores
getDoParWorkers()
setwd("~/Documents/Fluscape/serosolver/scripts/testing/coin_toss")
source("model_funcs.R")
source("probability_funcs.R")
source("proposal_funcs.R")
source("mcmc_funcs.R")
source("mcmc_funcs_marginal.R")

## Input parameters
n <- 20
coin_probs <- runif(n,0,0.2)
indivs <- 50
samps <- seq(1,n, by=1)
pars <- c(4, 0.3, 1)
fixed <- c(0,0,0)
fixed_probs <- rep(0,n)
#fixed_probs[1:2] <- 0
covMat_theta <- diag(length(fixed[which(fixed==0)]))
covMat_probs <- diag(length(fixed_probs[which(fixed_probs==0)]))
iter <- 50000

## Setup parameter names and simulated data
parNames <- c("boost","sigma","error")
coin_results <- sapply(coin_probs, function(x) sample(c(0,1),indivs,prob=c(1-x,x),replace=TRUE))
data_suggested_coins <- colSums(coin_results)/nrow(coin_results)
dat <- coin_toss_group(pars, coin_results)
dat <- measurement_error_group(pars,dat)

print(paste0("Coin probs: ", paste0(signif(coin_probs,3),collapse=" ")))
print(paste0("Data suggest coin values: ", paste0(data_suggested_coins,collapse=" ")))

startPars <- pars
startPars[1] <- runif(1,0,10)
startPars[2] <- runif(1,0,1)
startPars[3] <- runif(1,0,5)
startPars[which(fixed == 1)] <- pars[which(fixed == 1)]
startProbs <- runif(n,0,1)
startProbs[which(fixed_probs == 1)] <- coin_probs[which(fixed_probs == 1)]
start_coins <- matrix(sample(c(0,1),n*indivs,replace=TRUE),nrow=indivs)
res <- run_MCMC_group(startPars, startProbs, fixed, fixed_probs, coin_results,dat,samps, iter, 
                    covMat_theta, covMat_probs, thin=100,0.01,0.001,500,1,printF=1000,temp=1)
res1 <- run_MCMC_marginal(startPars, fixed, coin_results,dat,samps, iter, 
                      covMat_theta, thin=100,0.01,500,1,printF=1000,temp=1, samplePropn=0.1)

## Look at MCMC chain of process parameters
chain <- res[[2]]
#chain <- chain[,c(1,which(fixed == 0)+1)]
colnames(chain) <- c("sampno",parNames, paste0("coin.",1:n))
chain <- chain[chain[,"sampno"] > 0.1*iter,]
#plot(coda::as.mcmc(chain))


chain1 <- res1[[2]]
colnames(chain1) <- c("sampno",parNames)
chain1 <- chain1[chain1[,"sampno"] > 0.1*iter,]

#plot(coda::as.mcmc(chain1))
chains <- mcmc.list(as.mcmc(chain[,1:4]),as.mcmc(chain1))
plot(chains)
coin_chains <- mcmc.list(as.mcmc(res[[3]][,2:11]),as.mcmc(res1[[3]][,2:11]))
plot(coin_chains)

chain <- melt(as.data.frame(chain),id.vars="sampno")
colnames(chain) <- c("sampno","variable","value")
real_pars <- data.frame(variable=c(parNames,paste0("coin.",1:n)),y=c(pars,coin_probs),lower=c(3.5,rep(0,n+2)),upper=c(4.5,0.5,5,rep(1,n)))
fixed_pars <- c("sampno")
real_pars <- real_pars[!(real_pars$variable %in% fixed_pars),]
chain <- chain[!(chain$variable %in% fixed_pars),]

p_chain <- ggplot(chain) + geom_line(aes(x=sampno,y=value))+facet_wrap(~variable,scales="free_y")


## Look at posteriors vs. real parameters
p1 <- ggplot(chain) + 
  geom_density(aes(x=value,y=..scaled..)) +
  geom_vline(data=real_pars,aes(xintercept=y),col="red",linetype="dashed")+
  geom_blank(data=real_pars,aes(x=lower))+
  geom_blank(data=real_pars,aes(x=upper))+
  facet_wrap(~variable, scales="free_x")
 
## Coin flip outcome chain
coin_chain <- as.data.frame(res[[3]])
colnames(coin_chain) <- c("sampno",1:n,"indiv")
## First 10% burn in
coin_chain <- coin_chain[coin_chain$sampno > 0.1*iter,]
## Melt chain by sample number
melted_coin_chain <- melt(coin_chain[,1:(ncol(coin_chain)-1)], id.vars="sampno")
## Get proportion heads at each sample for each coin
melted_coin_chain <- ddply(melted_coin_chain, .(sampno, variable), function(x) sum(x$value)/indivs)
colnames(melted_coin_chain)[3] <- "value"
## Expand such that we have a column for each coin
#casted_coin_chain <- dcast(melted_coin_chain, sampno~variable)
## Get CI for each chain
casted_coin_chain <- ddply(melted_coin_chain, .(variable), function(x) quantile(x$value, c(0.025,0.5,0.975)))

tmp <- rbind(chain, melted_coin_chain)
tmp1 <- dcast(tmp, sampno~variable)

inferred_coins <- plyr::ddply(coin_chain, ~indiv, colMeans)
inferred_coins <- inferred_coins[,2:ncol(inferred_coins)]
mean_heads_indiv <- melt(inferred_coins,id.vars="indiv")
real_heads_indiv <- melt(coin_results)
colnames(real_heads_indiv) <- c("indiv","variable","value")
real_heads_indiv <- real_heads_indiv[real_heads_indiv$value == 1,]
p2 <- ggplot(mean_heads_indiv[mean_heads_indiv$indiv %in% 1:10,]) + 
  geom_point(aes(x=variable,y=value)) + 
  geom_vline(data=real_heads_indiv[real_heads_indiv$indiv %in% 1:10,],aes(xintercept=variable),col="red",linetype="dashed",alpha=0.4) +
  facet_wrap(~indiv)

## Actual coin weightings
real_coins <- data.frame(variable=1:n, y=coin_probs)
## MLE from actual coin flips
real_coins_data <- data.frame(variable=1:n, y=colSums(coin_results)/nrow(coin_results))
coin_prob_chain <- chain[chain$variable %in% paste0("coin.",1:n),]
colnames(coin_prob_chain) <- c("x","variable","value")
casted_coin_chain$variable <- as.integer(casted_coin_chain$variable)
coin_prob_chain$variable <- factor(coin_prob_chain$variable, levels=unique(as.character(coin_prob_chain$variable)))
coin_prob_chain$variable <- as.integer(coin_prob_chain$variable)
coin_prob_chain <- ddply(coin_prob_chain, ~variable, function(x) quantile(x$value, c(0.025,0.5,0.975)))

samp_tmp <- data.frame(xmin=samps-0.5,xmax=samps+0.5,ymin=0,ymax=1)
p3 <- ggplot() + 
  geom_rect(data=samp_tmp,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="gray80"),alpha=0.5)+
  geom_pointrange(data=coin_prob_chain,aes(x=variable-0.2,ymin=`2.5%`,y=`50%`,ymax=`97.5%`,col="coin prob chain")) + 
  geom_pointrange(data=casted_coin_chain,aes(x=variable+0.2,ymin=`2.5%`,y=`50%`,ymax=`97.5%`,col="coin chain")) + 
  geom_point(data=real_coins, aes(x=variable,y=y,col="real coin prob")) +
  geom_point(data=real_coins_data,aes(x=variable,y=y,col="real coin data")) +
  theme_bw() +
  scale_colour_manual(name="Key",values=c("coin prob chain"="blue","coin chain"="purple",
                                          "real coin prob"="black","real coin data"="red"),
                                  labels=c("Coin chain","Coin prob chain","Real coin data","Real coin prob")) +
  scale_fill_identity(name = 'Data', guide = 'legend',labels = c('Observed')) +
  xlab("Coin") +
  ylab("Inferred/known probability of heads")+
  scale_x_continuous(expand=c(0,0),limits=c(0.5,n+0.5),breaks=seq(1,n,by=1),labels=seq(1,n,by=1))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1))

p1
p2
p3
p_chain
