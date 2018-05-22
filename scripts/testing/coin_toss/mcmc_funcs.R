run_MCMC_group <- function(pars, probs,
                           fixed_pars, fixed_probs,
                           coin_results, dat, samps,
                           iter,covMat_theta, covMat_probs, 
                           thin=10,
                           step_theta,
                           step_prob,
                           adapt_freq,
                           adaptive_period,
                           printF=100,
                           temp=1){
  ## Some tuning parameters
  opt_tuning <- 0.2
  w <- 0.9
  
  ## Empty chains and storage objects
  n_indiv <- nrow(coin_results)
  liks <- numeric(iter/thin)
  chain <- matrix(nrow=iter/thin,ncol=length(pars) + length(probs) +1)
  coin_chain <- matrix(nrow=iter*n_indiv/thin,ncol=ncol(coin_results)+2)
  opt_chain <- matrix(nrow=adaptive_period,ncol=length(pars) + length(probs))
  
  proposed_coin_results <- coin_results
  probabs <- likelihood_group(pars, coin_results, dat, samps)/temp + hyper_prior_group(probs, coin_results)/temp
  probab <- sum(probabs) + prior(pars, probs)
  liks[1] <- probab
  
  
  chain[1,1] <- 1
  chain[1,2:(length(pars)+1)] <- pars
  chain[1, (length(pars)+2):ncol(chain)] <- probs
  
  opt_chain[1,] <- c(pars, probs)
  
  coin_chain[1:n_indiv,1] <- 1
  coin_chain[1:n_indiv,2:(ncol(coin_chain)-1)] <- coin_results
  coin_chain[1:n_indiv,ncol(coin_chain)] <- 1:n_indiv
  
  covMat0_theta <- covMat_theta
  covMat0_probs <- covMat_probs
  
  
  accepted_theta <- iter_theta <- accepted_theta_total <- iter_theta_total <-0
  accepted_probs <- iter_probs <- accepted_probs_total <- iter_probs_total <- 0
  accepted_coin <- iter_coin <- numeric(n_indiv)
  
  fixed <- which(fixed == 0)
  fixed_probs <- which(fixed_probs==0)
  
  
  
  index <- 1
  for(i in 2:iter){
    if(i %% 2 == 0){
      proposed <- proposal_theta(pars, fixed, covMat_theta*step_theta, covMat0_theta*step_theta)
      iter_theta <- iter_theta + 1
      iter_theta_total <- iter_theta_total + 1
      new_probabs <- likelihood_group(proposed, coin_results, dat, samps)/temp + hyper_prior_group(probs, coin_results)/temp
      new_probab <- sum(new_probabs) + prior(proposed, probs)
      #if(!is.finite(new_probab)) new_probab <- -100000000
      log_prob <- min(new_probab - probab,0)
      if(is.finite(log_prob) & log(runif(1))<log_prob){
        pars <- proposed
        probab <- new_probab
        probabs <- new_probabs
        accepted_theta <- accepted_theta + 1
        accepted_theta_total <- accepted_theta_total + 1
        #
      }
     # step_theta <- rm_scale(step_theta, i, log_prob, adaptive_period)
    }else if(i %% 3 == 0){
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
    if(i %% printF == 0){
      message(paste0("Iteration: ",i))
      message(paste0("Acceptance rate on theta: ", signif(accepted_theta/iter_theta,3)))
      message(paste0("Acceptance rate on probs: ", signif(accepted_probs/iter_probs,3)))
      message(paste0("Acceptance rate on coins: ", paste(signif(accepted_coin/iter_coin,3),collapse=" ")))
      message(paste0("Step size theta: ", step_theta))
      message(paste0("Step size prob: ", step_prob))
      message("\n")
      accepted_theta <- iter_theta <- 0
      accepted_probs <- iter_probs <- 0
      accepted_coin <- iter_coin <- numeric(n_indiv)
    }
    if(i %% thin == 0){
      liks[index] <- probab
      chain[index,1] <- i
      chain[index,2:(length(pars)+1)] <- pars
      chain[index, (length(pars)+2):ncol(chain)] <- probs
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),1] <- i
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),2:(ncol(coin_chain)-1)] <- coin_results
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),ncol(coin_chain)] <- 1:n_indiv
      index <- index + 1
    }
    
    if(i < adaptive_period){
      opt_chain[i,] <- c(pars, probs)
      if(i %% adapt_freq == 0){
        if(i > 0.1*adaptive_period){
          message("Adapting covariance matrix")
          message("\n")
          oldcovMat_theta <- covMat_theta
          oldcovMat_probs <- covMat_probs
          covMat_theta <- w*cov(opt_chain[1:i,fixed]) + (1-w)*oldcovMat_theta
          covMat_probs <- w*cov(opt_chain[1:i,(length(pars)+1):ncol(opt_chain)]) + (1-w)*oldcovMat_probs
          if(i > 0.9*adaptive_period){
            pcur_theta <- accepted_theta_total/iter_theta_total
            pcur_prob <- accepted_probs_total/iter_probs_total
            step_theta =max(0.00001,min(1,exp(log(step_theta)+(pcur_theta-0.234)*0.999^i)))
            step_prob =max(0.00001,min(1,exp(log(step_theta)+(pcur_prob-0.234)*0.999^i)))
            #step_theta <- scaletuning(step_theta, 0.234, pcur_theta)
            #step_prob <- scaletuning(step_prob, 0.234, pcur_prob)
          }
        } else {
          #pcur_theta <- accepted_theta_total/iter_theta_total
          #pcur_prob <- accepted_probs_total/iter_probs_total
          #step_theta =max(0.00001,min(1,exp(log(step_theta)+(pcur_theta-0.234)*0.999^i)))
          #step_theta <- scaletuning(step_theta, 0.234, pcur_theta)
          #step_prob <- scaletuning(step_prob, 0.234, pcur_prob)
      }
    }
   }
  }
  return(list(liks, chain, coin_chain))
}

rm_scale <- function(step_scale, mc, log_prob, adaptive_period)
{
  dd <- exp(log_prob)
  if( dd < -30 ){ dd <- 0 }
  dd <- min( dd, 1 )
  
  rm_temp <- ( dd - 0.23 )/( (mc+1)/(0.01*adaptive_period+1) )
  
  out <- step_scale*exp(rm_temp)
  
  out <- max( out, 0.02 )
  out <- min( out, 2)
  out
}
scaletuning <- function(step, popt,pcur){
  if(pcur ==1) pcur <- 0.99
  if(pcur == 0) pcur <- 0.01
  step = (step*qnorm(popt/2))/qnorm(pcur/2)
  if(step > 25) step <- 25
  step <- max(0.00001, step)
  return(step)
}
