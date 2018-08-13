run_MCMC_gibbs <- function(pars, 
                           fixed_pars,
                            coin_results, 
                            dat, 
                            samps,
                            iter,
                            covMat_theta,
                            thin=10,
                            step_theta,
                            lower_bounds=c(0,0,0),
                            upper_bounds=c(10,1,5),
                            adapt_freq,
                            adaptive_period,
                            printF=100,
                            temp=1,
                            sampPropn=0.5,
                           yearPropn=0.5,
                           swapPropn=0.5,
                           theta_proposal="univariate_theta",
                           Z_proposal = "gibbs",
                           alpha=1,
                           beta=1){
    opt_tuning <- 0.2
    w <- 0.9
    n_indiv <- nrow(coin_results)
  n_years <- ncol(coin_results)
  
  index <- 1
  liks <- numeric((iter)/thin +1)
  chain <- matrix(nrow=1 + (iter)/thin,ncol=length(pars)+1)
  coin_chain <- matrix(nrow=n_indiv + (iter*n_indiv)/thin,ncol=ncol(coin_results)+2)
  opt_chain <- matrix(nrow=n_indiv + adaptive_period,ncol=length(pars))
  
  coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),1] <- 1
  coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),2:(ncol(coin_chain)-1)] <- coin_results
  coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),ncol(coin_chain)] <- 1:n_indiv

  probabs <- likelihood_group(pars, coin_results, dat, samps)
  probab <- sum(probabs) + prior_gibbs(pars)  + inf_mat_prior(coin_results, alpha, beta)
  
  liks[1] <- probab
  chain[1,1] <- 1
  chain[1,2:(ncol(chain))] <- pars
  
  opt_chain[1,] <- pars
  
  accepted_theta_total <- iter_theta_total <-0
  if(theta_proposal == "univariate_theta"){
    accepted_theta <- iter_theta <- numeric(length(pars))
  } else {
    accepted_theta <- iter_theta <- accepted_theta_total <- iter_theta_total <-0
  }
  accepted_coin <- iter_coin <- 0
  fixed <- which(fixed_pars== 0)
  
  covMat0_theta <- covMat_theta
  
  ii <- 1
  index <- 2
  for(i in 2:iter){
    if(i %% 2 != 0){
      if(theta_proposal == "univariate_theta"){
        proposed <- univ_proposal(pars, lower_bounds, upper_bounds, step_theta, ii)
        iter_theta[ii] <- iter_theta[ii] + 1
      } else {
        proposed <- proposal_theta(pars, fixed, covMat_theta*step_theta, covMat0_theta*step_theta)
        iter_theta <- iter_theta + 1    
      }
      new_probabs <- likelihood_group(proposed, coin_results, dat, samps)
      new_probab <- sum(new_probabs) + prior_gibbs(proposed)  + inf_mat_prior(coin_results, alpha, beta)
      log_prob <- min(new_probab - probab,0)
      
      iter_theta_total <- iter_theta_total + 1
      
      if(is.finite(log_prob) & log(runif(1))<log_prob){
        pars <- proposed
        probab <- new_probab
        probabs <- new_probabs
        ## Which proposal on theta are we using?
        if(theta_proposal == "univariate_theta"){
          accepted_theta[ii] <- accepted_theta[ii] + 1
        } else {
          accepted_theta <- accepted_theta + 1 
        }
        accepted_theta_total <- accepted_theta_total + 1
      }
      ii <- ii + 1
      if(ii > length(pars)) ii <- 1
    } else {
      if(Z_proposal != "gibbs"){
        sampledI <- sample(1:n_indiv, floor(n_indiv*sampPropn))
        k <- floor(n_years*yearPropn)
        ## Choose which proposal on infection histories we're using
        if(Z_proposal == "simple"){
          proposed_coin_results <- coin_proposal_simple(coin_results,k)
        } else if(Z_proposal == "group"){
          proposed_coin_results <- coin_proposal_symmetric_group(coin_results,k,sampledI,swapPropn)
        } else {
          proposed_coin_results <- coin_proposal_by_year(coin_results,k,n_indiv,swapPropn)
        }
        new_probabs <- likelihood_group(pars, proposed_coin_results, dat, samps)
        new_probab <- sum(new_probabs) + prior_gibbs(pars)  + inf_mat_prior(proposed_coin_results, alpha, beta)
        log_prob <- min(new_probab - probab,0)
        iter_coin <- iter_coin + 1
        if(is.finite(log_prob) & log(runif(1))<log_prob){
          coin_results <- proposed_coin_results
          probab <- new_probab
          probabs <- new_probabs
          accepted_coin <- accepted_coin + 1
        }
      } else {
        coin_results <- coin_proposal_gibbs(pars, coin_results, dat, alpha, beta, indiv_propn=sampPropn, year_propn=yearPropn,swapPropn=swapPropn)
        iter_coin <- iter_coin + 1
        accepted_coin <- accepted_coin + 1
        probabs <- likelihood_group(pars, coin_results, dat, samps)
        probab <- sum(probabs) + prior_gibbs(pars) + inf_mat_prior(coin_results,alpha,beta)
      }
     
    }
    if(i %% printF == 0){
      message(cat("Iteration: ",i,sep="\t"))
      message(cat("Acceptance rate on theta: ", signif(accepted_theta/iter_theta,3),sep="\t"))
      message(cat("Acceptance rate on coins: ", signif(accepted_coin/iter_coin,3),sep="\t"))
      message(cat("Step size theta: ", step_theta,sep="\t"))
      message("\n")
    }
    
    if(i %% thin == 0){
      liks[index] <- probab
      chain[index,1] <- i
      chain[index,2:(length(pars)+1)] <- pars
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),1] <- i
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),2:(ncol(coin_chain)-1)] <- coin_results
      coin_chain[((index-1)*n_indiv + 1):(index*n_indiv),ncol(coin_chain)] <- 1:n_indiv
      index <- index + 1
    }
    if(i < adaptive_period){
      opt_chain[i,] <- pars
      if(i %% adapt_freq == 0){
        pcur_theta <- accepted_theta/iter_theta
        if(theta_proposal == "univariate_theta"){
          for(x in 1:length(step_theta)) step_theta[x] <- scaletuning1(step_theta[x], 0.44, pcur_theta[x])
        } else {
          ## Update covariance matrix
          if(i > 0.1*adaptive_period){
            message("Adapting covariance matrix")
            message("\n")
            oldcovMat_theta <- covMat_theta
            covMat_theta <- w*cov(opt_chain[1:i,fixed]) + (1-w)*oldcovMat_theta
            
            ## Update step size on covariance matrix
            if(i > 0.8*adaptive_period){
              step_theta <- scaletuning1(step_theta, 0.234, pcur_theta)
            }
          }
        }
        message(cat("Acceptance rate on theta: ", signif(accepted_theta/iter_theta,3),sep="\t"))
        message(cat("Acceptance rate on coins: ", paste(signif(accepted_coin/iter_coin,3),collapse=" "),sep="\t"))
        message(cat("Step size theta: ", step_theta,sep="\t"))
        message("\n")
        if(theta_proposal == "univariate_theta"){
          accepted_theta <- iter_theta <- numeric(length(pars))
        } else {  
          accepted_theta <- iter_theta <- 0
        }
        accepted_coin <- iter_coin <- 0
      }
    }
  }
  return(list(liks, chain, coin_chain))
}
