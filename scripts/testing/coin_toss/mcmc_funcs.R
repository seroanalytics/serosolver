run_MCMC_group <- function(pars, 
                           probs,
                           fixed_pars, 
                           fixed_probs,
                           coin_results, 
                           dat, 
                           samps,
                           iter,
                           lower_bounds=c(0,0,0),
                           upper_bounds=c(10,1,5),
                           covMat_theta, 
                           covMat_probs, 
                           thin=10,
                           step_theta,
                           step_prob,
                           adapt_freq,
                           adaptive_period,
                           printF=100,
                           temp=1,
                           sampPropn=0.5,
                           yearPropn=0.5,
                           alpha=1,
                           beta=1,
                           theta_proposal="univariate",
                           Z_proposal="simple"){
  ## Some tuning parameters
  opt_tuning <- 0.2
  w <- 0.9
  
  ## Empty chains and storage objects
  n_indiv <- nrow(coin_results)
  n_years <- ncol(coin_results)
  

  liks_all <- liks_data <- liks_prior <- liks <- numeric(1 + iter/thin)
  chain <- matrix(nrow=1 + iter/thin,ncol=length(pars) + length(probs) +1)
  coin_chain <- matrix(nrow=iter*n_indiv/thin,ncol=ncol(coin_results)+2)
  opt_chain <- matrix(nrow=n_indiv + adaptive_period,ncol=length(pars) + length(probs))
  
  proposed_coin_results <- coin_results
  probabs <- likelihood_group(pars, coin_results, dat, samps) + hyper_prior_group(probs, coin_results)
  probab <- sum(probabs) + prior_probs(pars, probs)
  liks[1] <- probab

  liks_data[1] <- sum(likelihood_group(pars, coin_results, dat, samps))
  liks_prior[1] <- sum(hyper_prior_group(probs, coin_results))
  liks_all[1] <- liks_data[1] + liks_prior[1] + prior_probs(pars, probs)
  
  chain[1,1] <- 1
  chain[1,2:(length(pars)+1)] <- pars
  chain[1, (length(pars)+2):ncol(chain)] <- probs
  
  opt_chain[1,] <- c(pars, probs)
  
  index <- 1
  coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),1] <- 1
  coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),2:(ncol(coin_chain)-1)] <- coin_results
  coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),ncol(coin_chain)] <- 1:n_indiv
  
  covMat0_theta <- covMat_theta
  covMat0_probs <- covMat_probs
  
  #accepted_theta <- iter_theta <- accepted_theta_total <- iter_theta_total <-0
  #accepted_probs <- iter_probs <- accepted_probs_total <- iter_probs_total <- 0
  accepted_theta_total <- iter_theta_total <-0
  accepted_probs_total <- iter_probs_total <- 0
  if(theta_proposal == "univariate"){
    accepted_theta <- iter_theta <- numeric(length(pars))
    accepted_probs <- iter_probs <- numeric(length(probs))
  } else {
    accepted_theta_total <- iter_theta_total <-0
    accepted_probs_total <- iter_probs_total <- 0
  }
  accepted_coin <- iter_coin <- numeric(n_indiv)
  
  fixed <- which(fixed == 0)
  fixed_probs <- which(fixed_probs==0)
  
  lower_bounds_probs <- rep(0, n_years)
  upper_bounds_probs <- rep(1,n_years)
  
  ii <- jj <- 1
  
  index <- 1
  for(i in 2:iter){
    if(i %% 2 == 0){
      #proposed <- proposal_theta(pars, fixed, covMat_theta*step_theta, covMat0_theta*step_theta)
      if(theta_proposal == "univariate"){
        proposed <- univ_proposal(pars, lower_bounds, upper_bounds, step_theta, ii)
        iter_theta[ii] <- iter_theta[ii] + 1
      } else {
        proposed <- proposal_theta(pars, fixed, covMat_theta*step_theta, covMat0_theta*step_theta)
        iter_theta <- iter_theta + 1    
      }
      iter_theta_total <- iter_theta_total + 1  
      new_probabs <- likelihood_group(proposed, coin_results, dat, samps) + hyper_prior_group(probs, coin_results)
      new_probab <- sum(new_probabs) + prior_probs(proposed, probs)
      log_prob <- min(new_probab - probab,0)
      if(is.finite(log_prob) & log(runif(1))<log_prob/temp){
        pars <- proposed
        probab <- new_probab
        probabs <- new_probabs
        
        ## Which proposal on theta are we using?
        if(theta_proposal == "univariate"){
          accepted_theta[ii] <- accepted_theta[ii] + 1
        } else {
         accepted_theta <- accepted_theta + 1 
        }
        accepted_theta_total <- accepted_theta_total + 1
      }
      ii <- ii + 1
      if(ii > length(pars)) ii <- 1
    }else if(i %% 3 == 0){
      ## Which proposal on probs are we using?
      if(theta_proposal == "univariate"){
        proposed_probs <- univ_proposal(probs, lower_bounds_probs, upper_bounds_probs, step_prob, jj)
        iter_probs[jj] <- iter_probs[jj] + 1
      } else {
        proposed_probs <- proposal_theta(probs, fixed_probs, covMat_probs*step_prob, covMat0_probs*step_prob)
        iter_probs <- iter_probs + 1
      }
      #proposed_coins <- coin_proposal_probs(coin_results, proposed_probs, 0.1*ncol(coin_results),0.1*nrow(coin_results))
      iter_probs_total <- iter_probs_total + 1
      new_probabs <- likelihood_group(pars, coin_results, dat, samps) + hyper_prior_group(proposed_probs, coin_results)
      new_probab <- sum(new_probabs) + prior_probs(pars, proposed_probs)

      log_prob <- min(new_probab - probab,0)
      if(is.finite(log_prob) & log(runif(1))<log_prob/temp){
        probs <- proposed_probs
        probab <- new_probab
        probabs <- new_probabs
        if(theta_proposal == "univariate"){
          accepted_probs[jj] <- accepted_probs[jj] + 1
        } else {
          accepted_probs <- accepted_probs + 1
        }
        accepted_probs_total <- accepted_probs_total + 1
      }
      jj <- jj + 1
      if(jj > length(probs)) jj <- 1
    } else {
      sampledI <- sample(1:n_indiv, floor(n_indiv*sampPropn))
      k <- floor(n_years*yearPropn)
      
      ## Choose which proposal on infection histories we're using
      if(Z_proposal == "simple"){
        proposed_coin_results <- coin_proposal_simple(coin_results,k,swapPropn)
      } else if(Z_proposal == "group"){
        proposed_coin_results <- coin_proposal_symmetric_group(coin_results,k,sampledI,swapPropn)
      } else if(Z_proposal == "years"){
        proposed_coin_results <- coin_proposal_by_year(coin_results,k,n_indiv,swapPropn)
      } else {
        proposed_coin_results <- coin_proposal_probs(coin_results, probs, k, sampledI,swapPropn)
      }
     
      new_probabs <- likelihood_group(pars, proposed_coin_results, dat, samps) + hyper_prior_group(probs, proposed_coin_results)
      iter_coin[sampledI] <- iter_coin[sampledI] + 1
      log_probs <- (new_probabs[sampledI] - probabs[sampledI])
      log_probs[log_probs > 0] <- 0
      x <- which(log(runif(length(sampledI))) < log_probs/temp)
      changeI <- sampledI[x]
      coin_results[changeI,] <- proposed_coin_results[changeI,]
      probabs[changeI] <- new_probabs[changeI]
      accepted_coin[changeI] <- accepted_coin[changeI] + 1
      probab <- sum(probabs) + prior_probs(pars, probs)
    }
    if(i %% printF == 0){
      message(cat("Iteration: ",i))
      message(cat("Acceptance rate on theta: ", signif(accepted_theta/iter_theta,3),sep="\t"))
      message(cat("Acceptance rate on probs: ", signif(accepted_probs/iter_probs,3),sep="\t"))
      message(cat("Acceptance rate on coins: ", paste(signif(accepted_coin/iter_coin,3),collapse=" "),sep="\t"))
      message("\n")
    
    }
    if(i %% thin == 0){
      liks[index] <- probab
      liks_data[index] <- sum(likelihood_group(pars, coin_results, dat, samps))
      liks_prior[index] <- sum(hyper_prior_group(probs, coin_results))
      liks_all[index] <- liks_data[index] + liks_prior[index] + prior_probs(pars, probs)
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
        pcur_theta <- accepted_theta/iter_theta
        pcur_prob <- accepted_probs/iter_probs
        if(theta_proposal == "univariate"){
          for(x in 1:length(step_theta)) step_theta[x] <- scaletuning1(step_theta[x], 0.44, pcur_theta[x])
          for(x in 1:length(step_prob)) step_prob[x] <- scaletuning1(step_prob[x], 0.44, pcur_prob[x])
        } else {
          ## Update covariance matrix
          if(i > 0.1*adaptive_period){
            message("Adapting covariance matrix")
            message("\n")
            oldcovMat_theta <- covMat_theta
            oldcovMat_probs <- covMat_probs
            covMat_theta <- w*cov(opt_chain[1:i,fixed]) + (1-w)*oldcovMat_theta
            covMat_probs <- w*cov(opt_chain[1:i,(length(pars)+1):ncol(opt_chain)]) + (1-w)*oldcovMat_probs
           
            ## Update step size on covariance matrix
            if(i > 0.8*adaptive_period){
              #step_theta =max(0.00001,min(1,exp(log(step_theta)+(pcur_theta-0.234)*0.999^i)))
              #step_prob =max(0.00001,min(1,exp(log(step_theta)+(pcur_prob-0.234)*0.999^i)))
              step_theta <- scaletuning1(step_theta, 0.234, pcur_theta)
              step_prob <- scaletuning1(step_prob, 0.234, pcur_prob)
            }
          }
        }
          message(cat("Acceptance rate on theta: ", signif(accepted_theta/iter_theta,3),sep="\t"))
          message(cat("Acceptance rate on probs: ", signif(accepted_probs/iter_probs,3),sep="\t"))
          message(cat("Acceptance rate on coins: ", paste(signif(accepted_coin/iter_coin,3),collapse=" "),sep="\t"))
          message(cat("Step size theta: ", step_theta,sep="\t"))
          message(cat("Step size prob: ", step_prob,sep="\t"))
          message("\n")
          if(theta_proposal == "univariate"){
            accepted_theta <- iter_theta <- numeric(length(pars))
            accepted_probs <- iter_probs <- numeric(length(probs))
          } else {  
            accepted_theta <- iter_theta <- 0
            accepted_probs <- iter_probs <- 0  
          }
          accepted_coin <- iter_coin <- numeric(n_indiv)
      }
    }
  }
  return(list(liks, chain, coin_chain,liks_data, liks_prior, liks_all))
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
