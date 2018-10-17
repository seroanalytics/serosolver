run_MCMC_marginal <- function(pars, 
                              fixed_pars,
                              coin_results, 
                              dat, 
                              samps,
                              iter,
                              covMat_theta,
                              thin=10,
                              step_theta,
                              adapt_freq,
                              adaptive_period,
                              printF=100,
                              temp=1,
                              samplePropn=0.5,
                              yearPropn=0.5,
                              alpha=1,
                              beta=1){
  ## Some tuning parameters
  opt_tuning <- 0.2
  w <- 0.9
  
  ## Empty chains and storage objects
  n_indiv <- nrow(coin_results)
  n_years <- ncol(coin_results)
  
  liks <- numeric(iter/thin)
  chain <- matrix(nrow=iter/thin,ncol=length(pars)+1)
  coin_chain <- matrix(nrow=iter*n_indiv/thin,ncol=ncol(coin_results)+2)
  opt_chain <- matrix(nrow=adaptive_period,ncol=length(pars))
  
  proposed_coin_results <- coin_results
  probabs <- likelihood_group(pars, coin_results, dat, samps)/temp
  probab <- sum(probabs) + prior(pars) + inf_mat_prior(coin_results, alpha, beta)
  liks[1] <- probab
  
  
  chain[1,1] <- 1
  chain[1,2:ncol(chain)] <- pars
  
  opt_chain[1,] <- c(pars)
  
  coin_chain[1:n_indiv,1] <- 1
  coin_chain[1:n_indiv,2:(ncol(coin_chain)-1)] <- coin_results
  coin_chain[1:n_indiv,ncol(coin_chain)] <- 1:n_indiv
  
  covMat0_theta <- covMat_theta
  
  accepted_theta <- iter_theta <- accepted_theta_total <- iter_theta_total <-0
  accepted_coin <- iter_coin <- 0
  
  fixed <- which(fixed == 0)
  
  
  index <- 1
  for(i in 2:iter){
    if(i %% 2 == 0){
      proposed <- proposal_theta(pars, fixed, covMat_theta*step_theta, covMat0_theta*step_theta)
      iter_theta <- iter_theta + 1
      iter_theta_total <- iter_theta_total + 1
      new_probabs <- likelihood_group(proposed, coin_results, dat, samps)/temp
      new_probab <- sum(new_probabs) + prior(proposed)  + inf_mat_prior(coin_results, alpha, beta)
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
    } else {
      sampledI <- sample(1:n_indiv, floor(n_indiv*samplePropn))
      k <- floor(n_years*yearPropn)
                                        # 1:n_indivmodl
      proposed_coin_results <- coin_proposal_symmetric_group(coin_results,k,sampledI)
      #proposed_coin_results <- coin_proposal_by_year(coin_results,1,floor(n_indiv*samplePropn))
      new_probabs <- likelihood_group(pars, proposed_coin_results, dat, samps)/temp
      new_probab <- sum(new_probabs) + prior(pars)  + inf_mat_prior(proposed_coin_results, alpha, beta)
      log_prob <- min(new_probab - probab,0)
      iter_coin <- iter_coin+ 1
      if(is.finite(log_prob) & log(runif(1))<log_prob){
        coin_results <- proposed_coin_results
        probab <- new_probab
        probabs <- new_probabs
        accepted_coin <- accepted_coin + 1
        #
      }
      
    }
    if(i %% printF == 0){
      message(paste0("Iteration: ",i))
      message(paste0("Acceptance rate on theta: ", signif(accepted_theta/iter_theta,3)))
      message(paste0("Acceptance rate on coins: ", signif(accepted_coin/iter_coin,3)))
      message(paste0("Step size theta: ", step_theta))
      message("\n")
      accepted_theta <- iter_theta <- 0
      accepted_coin <- iter_coin <- 0
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
    
    if(i < adaptive_period){
      opt_chain[i,] <- c(pars)
      if(i %% adapt_freq == 0){
        if(i > 0.1*adaptive_period){
          message("Adapting covariance matrix")
          message("\n")
          oldcovMat_theta <- covMat_theta
          covMat_theta <- w*cov(opt_chain[1:i,fixed]) + (1-w)*oldcovMat_theta
          if(i > 0.9*adaptive_period){
            pcur_theta <- accepted_theta_total/iter_theta_total
            step_theta =max(0.00001,min(1,exp(log(step_theta)+(pcur_theta-0.234)*0.999^i)))
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
