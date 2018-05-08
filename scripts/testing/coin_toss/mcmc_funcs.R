run_MCMC_group <- function(pars, fixed, dat, coin_results, iter,covMat, thin=10, samps){
  n_indiv <- nrow(coin_results)
  liks <- numeric(iter/thin)
  chain <- matrix(nrow=iter/thin,ncol=length(pars)+1)
  coin_chain <- matrix(nrow=iter*n_indiv/thin,ncol=ncol(coin_results)+2)
  
  proposed <- pars
  proposed_coin_results <- coin_results
  probabs <- likelihood_group(pars, coin_results, dat, samps) + hyper_prior_group(pars, coin_results)
  probab <- sum(probabs) + prior(pars)
  liks[1] <- probab
  chain[1,1] <- 1
  chain[1,2:ncol(chain)] <- pars
  coin_chain[1:n_indiv,1] <- 1
  coin_chain[1:n_indiv,2:(ncol(coin_chain)-1)] <- coin_results
  coin_chain[1:n_indiv,ncol(coin_chain)] <- 1:n_indiv
  accepted_theta <- iter_theta <- 0
  accepted_coin <- iter_coin <- numeric(n_indiv)
  fixed <- which(fixed == 0)
  index <- 1
  for(i in 2:iter){
    if(i %% thin == 0) message(i)
    if(i %% 2 == 0){
      proposed <- proposal_theta(pars, fixed, covMat, covMat)
      iter_theta <- iter_theta + 1
      new_probabs <- likelihood_group(proposed, coin_results, dat, samps) + hyper_prior_group(proposed, coin_results)
      new_probab <- sum(new_probabs) + prior(proposed)
      log_prob <- min(new_probab - probab,0)
      if(is.finite(log_prob) & log(runif(1))<log_prob){
        pars <- proposed
        probab <- new_probab
        probabs <- new_probabs
        accepted_theta <- accepted_theta + 1
      }
    } else {
      sampledI <- 1:n_indiv
      proposed_coin_results <- coin_proposal_symmetric_group(coin_results,1,sampledI)
      new_probabs <- likelihood_group(pars, proposed_coin_results, dat, samps) + hyper_prior_group(pars, proposed_coin_results)
      iter_coin[sampledI] <- iter_coin[sampledI] + 1
      log_probs <- new_probabs[sampledI] - probabs[sampledI]
      log_probs[log_probs > 0] <- 0
      x <- which(log(runif(length(sampledI))) < log_probs)
      changeI <- sampledI[x]
      coin_results[changeI,] <- proposed_coin_results[changeI,]
      probabs[changeI] <- new_probabs[changeI]
      accepted_coin[changeI] <- accepted_coin[changeI] + 1
      probab <- sum(probabs) + prior(pars)
    }
    if(i %% thin == 0){
      liks[index] <- probab
      chain[index,1] <- i
      chain[index,2:ncol(chain)] <- pars
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),1] <- i
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),2:(ncol(coin_chain)-1)] <- coin_results
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),ncol(coin_chain)] <- 1:n_indiv
      index <- index + 1
    }
  }
  message(paste0("Acceptance rate on theta: ", accepted_theta/iter_theta))
  message(paste0("Acceptance rate on coins: ", paste(accepted_coin/iter_coin,collapse=" ")))
  return(list(liks, chain, coin_chain))
}


run_MCMC_indiv <- function(pars, fixed, dat, coin_results, coin_probs, iter,covMat, thin=10){
  n_indiv <- nrow(coin_results)
  liks <- numeric(iter/thin)
  chain <- matrix(nrow=iter/thin,ncol=length(pars)+1)
  coin_prob_chain <- coin_chain <- matrix(nrow=iter*n_indiv/thin,ncol=ncol(coin_results)+2)
  proposed <- pars
  proposed_coin_results <- coin_results
  proposed_coin_probs <- coin_probs
  
  probabs <- likelihood_group(pars, coin_results, dat) + prior_group(coin_probs, coin_results)
  probab <- sum(probabs) + prior(pars)
  
  
  liks[1] <- probab
  chain[1,1] <- 1
  chain[1,2:ncol(chain)] <- pars
  coin_chain[1:n_indiv,1] <- 1
  coin_chain[1:n_indiv,2:(ncol(coin_chain)-1)] <- coin_results
  coin_chain[1:n_indiv,ncol(coin_chain)] <- 1:n_indiv
  accepted_theta <- iter_theta <- 0
  accepted_coin <- iter_coin <- numeric(n_indiv)
  fixed <- which(fixed == 0)
  index <- 1
  for(i in 2:iter){
    if(i %% thin == 0) message(i)
    if(i %% 3 == 0){
      proposed <- proposal_theta(pars, fixed, covMat, covMat)
      iter_theta <- iter_theta + 1
      new_probabs <- likelihood_group(proposed, coin_results, dat) + prior_group(coin_probs, coin_results)
      new_probab <- sum(new_probabs) + prior(proposed)
      log_prob <- min(new_probab - probab,0)
      if(is.finite(log_prob) & log(runif(1))<log_prob){
        pars <- proposed
        probab <- new_probab
        probabs <- new_probabs
        accepted_theta <- accepted_theta + 1
      }
    } else if(i %% 3 == 1){
      sampledI <- 1:n_indiv
      proposed_coin_probs <- coin_prob_proposal(coin_probs,1,sampledI)
      new_probabs <- likelihood_group(pars, coin_results, dat) + prior_group(proposed_coin_probs, coin_results)
      #iter_coin[sampledI] <- iter_coin[sampledI] + 1
      log_probs <- new_probabs[sampledI] - probabs[sampledI]
      log_probs[log_probs > 0] <- 0
      x <- which(log(runif(length(sampledI))) < log_probs)
      changeI <- sampledI[x]
      coin_probs[changeI,] <- proposed_coin_probs[changeI,]
      probabs[changeI] <- new_probabs[changeI]
      probab <- sum(probabs) + prior(pars)
    } else {
      sampledI <- 1:n_indiv
      proposed_coin_results <- coin_proposal_symmetric_group(coin_results,1,sampledI)
      new_probabs <- likelihood_group(pars, proposed_coin_results, dat) + prior_group(coin_probs, proposed_coin_results)
      iter_coin[sampledI] <- iter_coin[sampledI] + 1
      log_probs <- new_probabs[sampledI] - probabs[sampledI]
      log_probs[log_probs > 0] <- 0
      x <- which(log(runif(length(sampledI))) < log_probs)
      changeI <- sampledI[x]
      coin_results[changeI,] <- proposed_coin_results[changeI,]
      accepted_coin[changeI] <- accepted_coin[changeI] + 1
      probabs[changeI] <- new_probabs[changeI]
      probab <- sum(probabs) + prior(pars)
    }
    if(i %% thin == 0){
      liks[index] <- probab
      chain[index,1] <- i
      chain[index,2:ncol(chain)] <- pars
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),1] <- i
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),2:(ncol(coin_chain)-1)] <- coin_results
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),ncol(coin_chain)] <- 1:n_indiv
      
      coin_prob_chain[((index-1)*n_indiv + 1):(index*n_indiv),1] <- i
      coin_prob_chain[((index-1)*n_indiv + 1):(index*n_indiv),2:(ncol(coin_prob_chain)-1)] <- coin_probs
      coin_prob_chain[((index-1)*n_indiv + 1):(index*n_indiv),ncol(coin_prob_chain)] <- 1:n_indiv
      
      index <- index + 1
    }
  }
  message(paste0("Acceptance rate on theta: ", accepted_theta/iter_theta))
  message(paste0("Acceptance rate on coins: ", paste(accepted_coin/iter_coin,collapse=" ")))
  return(list(liks, chain, coin_chain, coin_prob_chain))
}