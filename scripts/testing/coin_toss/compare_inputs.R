## Starting conditions and run MCMC
kmax <- 5
tmps <- seq(0.1,1,length.out=kmax)
foreach(i=1:kmax) %dopar% {
  startPars <- runif(3,0,2)
  startPars[which(fixed == 1)] <- pars[which(fixed == 1)]
  startProbs <- runif(n,0,1)
  start_coins <- matrix(sample(c(0,1),n*indivs,replace=TRUE),nrow=indivs)
  res <- run_MCMC_group(startPars, startProbs, fixed, fixed_probs, coin_results,dat,samps, iter, 
                        covMat_theta, covMat_probs, thin=100,0.001,0.00005,500,1,printF=1000,tmps[i])
  write.csv(res[[1]],paste0("likelihood_",i,".csv"))
  write.csv(res[[2]],paste0("theta_",i,".csv"))
  write.csv(res[[3]],paste0("coins_",i,".csv"))
}
chains <- NULL
for(i in 1:kmax){
  chain <- read.csv(paste0("theta_",i,".csv"))
  chain <- as.mcmc(chain[chain[,"V1"] > 10000,3:ncol(chain)])
  chains[[i]] <- chain
}
chains <- as.mcmc.list(chains)
ess_no_adaptation <- effectiveSize(chains)

foreach(i=1:kmax) %dopar% {
  startPars <- runif(3,0,2)
  startPars[which(fixed == 1)] <- pars[which(fixed == 1)]
  startProbs <- runif(n,0,1)
  start_coins <- matrix(sample(c(0,1),n*indivs,replace=TRUE),nrow=indivs)
  res <- run_MCMC_group(startPars, startProbs, fixed, fixed_probs, start_coins,dat,samps, iter, 
                        covMat_theta, covMat_probs, thin=100,1,1,500,100000,printF=1000)
  write.csv(res[[1]],paste0("likelihood_",i,".csv"))
  write.csv(res[[2]],paste0("theta_",i,".csv"))
  write.csv(res[[3]],paste0("coins_",i,".csv"))
}
chains <- NULL
for(i in 1:kmax){
  chain <- read.csv(paste0("theta_",i,".csv"))
  chain <- as.mcmc(chain[chain[,"V1"] > 50000,3:ncol(chain)])
  chains[[i]] <- chain
}
chains <- as.mcmc.list(chains)
ess_adaptation <- effectiveSize(chains)


foreach(i=1:kmax) %dopar% {
  startPars <- runif(3,0,2)
  startPars[which(fixed == 1)] <- pars[which(fixed == 1)]
  startProbs <- runif(n,0,1)
  start_coins <- matrix(sample(c(0,1),n*indivs,replace=TRUE),nrow=indivs)
  res <- run_MCMC_group(startPars, startProbs, fixed, fixed_probs, start_coins,dat,samps, iter, 
                        covMat_theta, covMat_probs, thin=100,0.001,0.005,500,100000,printF=1000)
  write.csv(res[[1]],paste0("likelihood_",i,".csv"))
  write.csv(res[[2]],paste0("theta_",i,".csv"))
  write.csv(res[[3]],paste0("coins_",i,".csv"))
}
chains <- NULL
for(i in 1:kmax){
  chain <- read.csv(paste0("theta_",i,".csv"))
  chain <- as.mcmc(chain[chain[,"V1"] > 50000,3:ncol(chain)])
  chains[[i]] <- chain
}
chains <- as.mcmc.list(chains)
ess_adaptation_b <- effectiveSize(chains)

beepr::beep(4)


plot(chains)
