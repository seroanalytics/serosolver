else if(i %% 3 == 0){
  proposed_probs <- proposal_theta(probs, fixed_probs, covMat_probs*step_prob, covMat0_probs*step_prob)
  iter_probs <- iter_probs + 1
  iter_probs_total <- iter_probs_total + 1
  new_probabs <- likelihood_group(pars, coin_results, dat, samps)/temp + hyper_prior_group(proposed_probs, coin_results)/temp
  new_probab <- sum(new_probabs) + prior(pars, proposed_probs)
  #if(!is.finite(new_probab)) new_probab <- -100000000
  log_prob <- min(new_probab - probab,0)
  if(is.finite(log_prob) & log(runif(1))<log_prob){
    probs <- proposed_probs
    probab <- new_probab
    probabs <- new_probabs
    accepted_probs <- accepted_probs + 1
    accepted_probs_total <- accepted_probs_total + 1
    #
  }
  
  #step_prob <- rm_scale(step_prob, i, log_prob, adaptive_period)
} else {
  sampledI <- 1:n_indiv
  proposed_coin_results <- coin_proposal_symmetric_group(coin_results,1,sampledI)
  new_probabs <- likelihood_group(pars, proposed_coin_results, dat, samps)/temp + hyper_prior_group(probs, proposed_coin_results)
  iter_coin[sampledI] <- iter_coin[sampledI] + 1
  log_probs <- (new_probabs[sampledI] - probabs[sampledI])
  log_probs[log_probs > 0] <- 0
  x <- which(log(runif(length(sampledI))) < log_probs/temp)
  changeI <- sampledI[x]
  coin_results[changeI,] <- proposed_coin_results[changeI,]
  probabs[changeI] <- new_probabs[changeI]
  accepted_coin[changeI] <- accepted_coin[changeI] + 1
  probab <- sum(probabs) + prior(pars, probs)
}

n <- 3
coin_probs <- runif(n,0,0.2)
indivs <- 200
samps <- seq(1,n, by=1)
coin_results <- sapply(coin_probs, function(x) sample(c(0,1),indivs,prob=c(1-x,x),replace=TRUE))
iter <- 10000
tmp <- matrix(nrow=iter,ncol=ncol(coin_results))

tmp_coin_results <- coin_results
lambdas <- runif(iter)
p_lambda <- numeric(iter)
N <- indivs
alpha <- beta <- 1
for(i in 1:iter){
  #tmp_coin_results <- coin_proposal_probs(tmp_coin_results, coin_probs, n, indivs)
  tmp_coin_results <- coin_proposal_gibbs(tmp_coin_results, 1000, 1000)
  tmp[i,] <- colSums(tmp_coin_results)
  x <- sum(tmp_coin_results[,1])
  lambdas[i] <- rbeta(1, alpha + x, beta + N - x)
  #p_lambda[i] <- prod((lambdas[i])^tmp_coin_results[,1] *(1-lambdas[i])^tmp_coin_results[,1])*dbeta(lambdas[i],1,1)/(beta(x + alpha, N - x + beta)/beta(alpha, beta))
    #prod((lambdas[i])^tmp_coin_results[,1] *(1-lambdas[i])^tmp_coin_results[,1])*dbeta(lambdas[i],1,1)/(beta(x + alpha, N - x + beta)/beta(alpha, beta))
    
    #/exp(inf_mat_prior(tmp_coin_results, 1,1))
}  
#plot(density(rbinom(10000,indivs,0.5)))
plot(density(tmp[,1]/indivs),col="red")
plot(density(lambdas))
#plot(p_lambda~lambdas)
#plot(as.mcmc(tmp/indivs))
summary(as.mcmc(tmp[,1]/indivs))[[1]]
summary(as.mcmc(lambdas))[[1]]


library(data.table)
library(plyr)
chain <- read.csv("real_lambda_1_1_chain.csv")
ages <- read.csv("real_lambda_1_1_ages.csv")

chain <- chain[chain$sampno > 100000,]
infChain <- data.table::fread("real_lambda_1_1_infectionHistories.csv")
infChain <- infChain[infChain$sampno > 100000,]
setkey(infChain, "sampno","j")
infChain1 <- infChain[,list(V1=sum(x)),by=key(infChain)]

plot(density(chain$lambda.2))
lines(density(infChain1[infChain1$j==3,]$V1/n_alive[3]),col="red")



plot(density(chain$lambda.5))
lines(density(infChain1[infChain1$j==6,]$V1/n_alive[6]),col="red")
