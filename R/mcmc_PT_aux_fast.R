create_run_MCMC_single_iter_fn_fast<- function(unfixed_pars,unfixed_par_length,
                                               lower_bounds,upper_bounds,
                                               age_mask, strain_mask,
                                               posterior_simp, extra_probabilities, proposal_gibbs,
                                               alpha, beta, hist_sample_prob, n_infs_vec, swap_propn, move_size,
                                               n_alive, switch_sample, hist_switch_prob, year_swap_propn
                                               ){
    f <- function(i, par_i,
                  current_pars,infection_histories,
                  likelihoods, total_likelihood,
                  prior_prob, posterior,
                  tempaccepted,tempiter,
                  steps,temp){
        ## Whether to swap entire year contents or not - only applies to gibbs sampling
        inf_swap_prob <- runif(1)

        if(i %% switch_sample != 0){
            ## If using univariate proposals
            ## For each parameter (Gibbs)
            j <- unfixed_pars[par_i]
            par_i <- par_i + 1
            if(par_i > unfixed_par_length) par_i <- 1
            proposal <- univ_proposal(current_pars, lower_bounds, upper_bounds, steps,j)
            tempiter[j] <- tempiter[j] + 1
            ## If using multivariate proposals
            ## Propose new parameters and calculate posterior
            ## Check that all proposed parameters are in allowable range
            if(!any(
                    proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                    proposal[unfixed_pars] > upper_bounds[unfixed_pars]
                )
               ){
                ## Calculate new likelihood and find difference to old likelihood
                new_likelihoods <- posterior_simp(proposal,infection_histories)
                new_total_likelihood <- sum(new_likelihoods) # Total
                new_prior_prob <- extra_probabilities(proposal, infection_histories) # Prior
                new_posterior <- new_total_likelihood + new_prior_prob # Posterior

                log_prob <- new_posterior-posterior
                if(!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)){
                    log_prob <- min(log_prob, 0)                    
                    if(log(runif(1)) < log_prob/temp){
                        ## Accept with probability 1 if better, or proportional to
                        ## difference if not
                        current_pars <- proposal
                        ## Store acceptances
                        tempaccepted[j] <- tempaccepted[j] + 1
                        likelihoods <- new_likelihoods
                        prior_prob <- new_prior_prob
                        total_likelihood <- new_total_likelihood
                        posterior <- new_posterior
                    }
                }
            }
        } else {
            ## Need to temporarily store current parameters as new pars, as
            ## might change with lambda swap step
            n_indiv <- nrow(infection_histories)
            indivSubSample <- sample(1:n_indiv, ceiling(hist_sample_prob*n_indiv))
            proposal <- current_pars
            if(inf_swap_prob > hist_switch_prob){
                prop_gibbs <- proposal_gibbs(proposal, infection_histories,
                                             likelihoods,
                                             indivSubSample, 
                                             alpha, beta,
                                             n_infs_vec,swap_propn,move_size,
                                             temp)
                new_likelihoods <- prop_gibbs$old_probs
                newInfectionHistories <- prop_gibbs$new_infection_history
            } else {
                newInfectionHistories <- inf_hist_swap(infection_histories, age_mask, strain_mask,
                                                       year_swap_propn, move_size)
                new_likelihoods <- posterior_simp(proposal, newInfectionHistories)/temp
            }
            
            ## Calculate new likelihood with these infection histories
            new_total_likelihood <- sum(new_likelihoods)
            new_prior_prob <- extra_probabilities(proposal, newInfectionHistories)
            new_posterior <- new_total_likelihood + new_prior_prob
            log_prob <- new_posterior - posterior
            
            if(inf_swap_prob > hist_switch_prob){
                if(!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)){
                    infection_histories <- newInfectionHistories
                    likelihoods <- new_likelihoods
                    total_likelihood <- new_total_likelihood
                    prior_prob <- new_prior_prob
                    posterior <- new_posterior
                }
            } else {
                if(!identical(newInfectionHistories, infection_histories)){
                    log_prob <- new_posterior - posterior
                    if(!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)){
                        log_prob <- min(log_prob, 0)
                        if(log(runif(1)) < log_prob){
                            infection_histories <- newInfectionHistories
                            current_pars <- proposal
                            likelihoods <- new_likelihoods
                            total_likelihood <- new_total_likelihood
                            prior_prob <- new_prior_prob
                            posterior <- new_posterior                        
                        }
                    }
                }
            }            
            list("i"=i, "par_i" = par_i,
                 "current_pars" = current_pars,"infection_histories"=infection_histories,
                 "likelihoods"=likelihoods,"total_likelihood" = total_likelihood,
                 "prior_prob"=prior_prob,"posterior"=posterior,
                 "tempaccepted" = tempaccepted,"tempiter" = tempiter,
                 "steps"=steps,"temp"=temp)
        }
        f
    }


#' performs parallel tempering - Ada Yan
#' 
#' @param mcmc_list a list of lists: values, log likelihood etc. of parallel MCMC chains
#' @param temperatures numeric vector: temperatures of chains
#' @param offset integer: 0 or 1. 0 = swap chains 1 <-> 2, 3 <-> 4...
#' 1 = swap chains 2<->3, 4<->5...
#' 
#' @return a list of lists: values, log likelihood etc. of paralle chains after parallel tempering
#' @export
parallel_tempering <- function(mcmc_list, temperatures, offset){
    recorded_swaps <- double(length(mcmc_list) - 1)
    
    ## extract current probabilities and log likelihoods
    all_likelihood <- vapply(mcmc_list, function(x) x$total_likelihood, double(1))
    all_prior_prob<- vapply(mcmc_list, function(x) x$prior_prob, double(1))
    all_posterior <- vapply(mcmc_list, function(x) x$posterior, double(1))
    all_likelihoods <- lapply(mcmc_list, function(x) x$likelihoods)
    
    all_current_pars <- lapply(mcmc_list, function(x) x$current_pars)
    all_infection_histories <- lapply(mcmc_list, function(x) x$infection_histories)
    
    ## decide which chains to swap
    if((offset + 1) <= (length(mcmc_list) - 1)){
        swap_ind <- seq(offset + 1, length(mcmc_list) - 1, by = 2)
        
        decide_if_swap <- function(x,y){
            delta <- (1 / temperatures[y] - 1 / temperatures[x]) *
                (all_posterior[x] - all_posterior[y])
            runif(1) <= exp(delta)
        }

        swaps <- vapply(swap_ind, function(x) decide_if_swap(x,x+1), logical(1))
        swap_ind <- swap_ind[swaps]

        ## perform swap
        perform_swap <- function(vec, swap_ind){
            vec_new <- vec
            vec_new[swap_ind] <- vec[swap_ind + 1]
            vec_new[swap_ind + 1] <- vec[swap_ind]
            vec_new
        }
        all_current_pars <- perform_swap(all_current_pars, swap_ind)
        all_likelihood <- perform_swap(all_likelihood, swap_ind)
        all_likelihoods <- perform_swap(all_likelihoods, swap_ind)
        all_prior_prob <- perform_swap(all_prior_prob, swap_ind)
        all_posterior <- perform_swap(all_posterior, swap_ind)
        all_infection_histories <- perform_swap(all_infection_histories, swap_ind)
        new_list <- Map(function(x,y,z,b,c,q) list(current_pars = x, infection_histories=y,
                                                   likelihoods=z, total_likelihood = b, 
                                                   prior_prob=c,posterior=q
                                                   ),
                        all_current_pars, all_infection_histories,
                        all_likelihoods, all_likelihood, 
                        all_prior_prob, all_posterior
                        )
        
        mcmc_list <- Map(modifyList, mcmc_list, new_list)
        recorded_swaps[swap_ind] <- 1
    }
    list("swaps" = recorded_swaps, "mcmc_list" = mcmc_list)
}

#' calibrate temperatures for parallel chains - Ada Yan
#'
#' @param temperatures vector of length n: current temperatures of chains
#' @param swap_ratio vector of length n - 1: (proportion of accepted swaps
#' out of proposed swaps)/2
#' the factor of 2 arises because of the way the swaps are recorded, and how
#' we alternate between swapping 1<->2, 3<->4.... and 2<->3, 4<->5...
#' @return vector of length n: new temperatures of chains
#'
calibrate_temperatures <- function(temperatures,swap_ratio) {
    diff_temp <- diff(temperatures)
    ## find chains between which the swap ratio is too large
    too_large = swap_ratio > .2 # note factor of 2 from main text -- see above
    ## find chains between which the swap ratio is too small
    too_small = swap_ratio < .05
    ## adjust differences between temperatures accordingly
    diff_temp = diff_temp*(too_large*1.5 + too_small*.75 + (1-too_large-too_small))
    ## reconstruct temperatures from their differences
    cumsum(c(temperatures[1],diff_temp))    
}
