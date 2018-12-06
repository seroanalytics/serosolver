    
    

run_MCMC_pt_fast <- function(par_tab,
                             titre_dat,
                             antigenic_map,
                             mcmcPars=c(),
                             startInfHist=NULL,
                             filename="test",
                             CREATE_POSTERIOR_FUNC=create_posterior_func_fast,
                             CREATE_PRIOR_FUNC=NULL,
                             seed,
                             mu_indices=NULL,
                             measurement_indices=NULL,
                             measurement_random_effects=FALSE,
                             OPT_TUNING=0.1,
                             solve_likelihood=TRUE,
                             n_alive=NULL,
                             ...){
    version <- 1
    
    ## Extract MCMC parameters
    mcmcPars_used <- list("iterations"=50000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,
                       "adaptive_period"=10000,
                       "save_block"=100,"thin2"=10,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
                       "inf_propn"=1, "moveSize"=5,"histOpt"=1,"swapPropn"=0.5,
                       "hist_switch_prob"=0,"year_swap_propn"=1,
                       "temperature"=seq(1,10,by=1),"parallel_tempering_iter"=5, "n_iter"=5)
    mcmcPars_used[names(mcmcPars)] <- mcmcPars

    ## Extract MCMC parameters
    iterations <- mcmcPars_used[["iterations"]] # How many iterations to run after adaptive period
    popt <- mcmcPars_used[["popt"]] # Desired optimal acceptance rate
    popt_hist <- mcmcPars_used[["popt_hist"]]
    opt_freq<- mcmcPars_used[["opt_freq"]] # How often to adjust step size
    thin <- mcmcPars_used[["thin"]] # Save only every nth iterations for theta sampling
    adaptive_period<- mcmcPars_used[["adaptive_period"]] # How many iterations for adaptive period
    save_block <- mcmcPars_used[["save_block"]] # How many post-thinning iterations to store before saving to disk
    histTabThin <- mcmcPars_used[["thin2"]] # Save only every nth iterations for infection history sampling
    histSampleProb <- mcmcPars_used[["histSampleProb"]] # What proportion of infection histories to sample each step
    switch_sample <- mcmcPars_used[["switch_sample"]] # Resample infection histories every n iterations
    burnin <- mcmcPars_used[["burnin"]] # Run this many iterations before attempting adaptation. Idea is to reduce getting stuck in local maxima
    moveSize <- mcmcPars_used[["moveSize"]] # Number of infections to move/remove/add in each proposal step
    inf_propn <- mcmcPars_used[["inf_propn"]] # Number of infections to move/remove/add in each proposal step
    nInfs <- floor(length(antigenic_map$inf_years)*inf_propn)
    histOpt <- mcmcPars_used[["histOpt"]] # Should infection history proposal step be adaptive?
    swapPropn <- mcmcPars_used[["swapPropn"]] # If using gibbs, what proportion of proposals should be swap steps?
    hist_switch_prob <- mcmcPars_used[["hist_switch_prob"]] # If using gibbs, what proportion of iterations should be swapping contents of two time periods?
    year_swap_propn <- mcmcPars_used[["year_swap_propn"]] # If gibbs and swapping contents, what proportion of these time periods should be swapped?
    temperatures <- mcmcPars[["temperature"]]
    parallel_tempering_iter <- mcmcPars[["parallel_tempering_iter"]]
    n_iter <- mcmcPars[["n_iter"]]
###################################################################

    ###############
    ## Currently only doing this for gibbs version
    ###############
    histProposal <- 2
    prior_on_total <- FALSE
    
    ## Extract current pars from first par_tab
    current_pars <- par_tab[[1]]$values

    ## Starting pars and step sizes for each rung seperately -
    ## this is the only bit we need for each rung
    start_pars <- lapply(par_tab, function(x) x$values)
    steps <- lapply(par_tab, function(x) x$steps)
    par_tab <- par_tab[[1]]

    ## All these bits the same for each rung
    param_length <- nrow(par_tab)
    unfixed_pars <- which(par_tab$fixed == 0)
    unfixed_par_length <- nrow(par_tab[par_tab$fixed== 0,])
    par_names <- as.character(par_tab$names)
    
    ## Parameter constraints
    lower_bounds <- par_tab$lower_bound
    upper_bounds <- par_tab$upper_bound
    fixed <- par_tab$fixed

    ## CURRENTLY NOT USED, AS USING GIBBS VERSION
    ## If using lambda terms, pull their indices out of the parameter table
    lambda_indices <- NULL
    if("lambda" %in% par_names){
        lambda_indices <- which(par_tab$names == "lambda")
    }
    alpha <- par_tab[par_tab$names == "alpha","values"]
    beta <- par_tab[par_tab$names == "beta","values"]

    ## Only doing univariate proposals
    tempaccepted <- tempiter <- integer(param_length)
    reset <- integer(param_length)
    reset[] <- 0

    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste0(filename,"_chain.csv")
    infectionHistory_file <- paste0(filename,"_infection_histories.csv")

    strain_isolation_times <- unique(antigenic_map$inf_years) # How many strains are we testing against and what time did they circulate
    samplingTimes <- unique(titre_dat$sample) # What are the range of sampling times?
    n_strain <- length(strain_isolation_times) # How many strains could an individual see?
    n_indiv <- length(unique(titre_dat$individual)) # How many individuals in the titre_dat?
    n_groups <- length(unique(titre_dat$group)) # How many groups in the titre_dat?
    individuals <- 1:n_indiv # Create vector of individuals
    groups <- 1:n_groups # Create vector of groups


 ###################
    ## Housekeeping for infection history chain
    ###################
    ## To store acceptance rate of entire time period infection history swaps
    infection_historySwapN <- infection_historySwapAccept <- 0

    histiter <- histaccepted <- histiter_add <- histaccepted_add <- histiter_move <- histaccepted_move <- histreset <- integer(n_indiv)
    
    nInfs_vec <- rep(nInfs, n_indiv) # How many infection history moves to make with each proposal
    moveSizes <- rep(moveSize, n_indiv) # How many years to move in smart proposal step

    
###############
    ## Create age mask
    ## -----------------------
    ## Note that ages for all groups must be from same reference point
    ## -----------------------
###############
    if(!is.null(titre_dat$DOB)){
        DOBs <- unique(titre_dat[,c("individual","DOB")])[,2]
    } else {
        DOBs <- rep(min(strain_isolation_times), n_indiv)
    }
    age_mask <- create_age_mask(DOBs, strain_isolation_times)
    ## Create strain mask
    strain_mask <- create_strain_mask(titre_dat,strain_isolation_times)
    masks <- data.frame(cbind(age_mask, strain_mask))
    ## Number of people that were born before each year and have had a sample taken since that year happened
    if(is.null(n_alive)) n_alive <- sapply(seq(1,length(strain_isolation_times)), function(x) nrow(masks[masks$age_mask <=x & masks$strain_mask >= x,]))   
    
    ## Create posterior calculating function
    posterior_simp <- protect(CREATE_POSTERIOR_FUNC(par_tab,
                                            titre_dat,
                                            antigenic_map,
                                            version,
                                            solve_likelihood,
                                            age_mask,
                                            measurement_indices_by_time=measurement_indices,
                                            mu_indices=mu_indices,
                                            n_alive=n_alive,
                                            function_type=1,
                                            ...))
    if(!is.null(CREATE_PRIOR_FUNC)){
        prior_func <- CREATE_PRIOR_FUNC(par_tab)
    }
    ## If using gibbs proposal on infection_history, create here
    proposal_gibbs <- CREATE_POSTERIOR_FUNC(par_tab,
                                            titre_dat,
                                            antigenic_map,
                                            version,
                                            solve_likelihood,
                                            age_mask,
                                            measurement_indices_by_time=measurement_indices,
                                            mu_indices=mu_indices,
                                            n_alive=n_alive,
                                            function_type=2,
                                            ...)

    ## If using random effects on mu, need to include hyperprior term on mu
    ## We can't do this in the main posterior function, because this term
    ## applies to the overall posterior whereas the main posterior function
    ## returns each individual's posterior
    if (!is.null(mu_indices)) {
        prior_mu <- create_prior_mu(par_tab)
    }
    if (measurement_random_effects) {
        prior_shifts <- create_prob_shifts(par_tab)
    }
######################
    
    ## Setup initial conditions
    ## I think this needs to be a list
    infection_histories = startInfHist[[1]]

    ## Initial likelihood
    likelihoods <- posterior_simp(current_pars,infection_histories)
    ## Initial total likelihood
    total_likelihood <- sum(likelihoods)
    n_alive_tot <- sum(n_alive)
    ## Create closure to add extra prior probabilities, to avoid re-typing later
    extra_probabilities <- function(prior_pars, prior_infection_history){
        prior_probab <- 0
        if(histProposal == 2){
            if(prior_on_total){
                prior_probab <- prior_probab + inf_mat_prior_total_cpp(prior_infection_history, n_alive_tot, alpha, beta)
            } else {
                prior_probab <- prior_probab + inf_mat_prior_cpp(prior_infection_history, n_alive, alpha, beta)
            }
        }
        if(!is.null(CREATE_PRIOR_FUNC)) prior_probab <- prior_probab + prior_func(prior_pars)
        if(!is.null(mu_indices)) prior_probab <- prior_probab + prior_mu(prior_pars)
        if(measurement_random_effects) prior_probab <- prior_probab + prior_shifts(prior_pars)
        prior_probab
    }
    ## Initial total prior prob
    prior_prob <- extra_probabilities(current_pars, infection_histories)

    ## Initial posterior prob
    posterior <- total_likelihood + prior_prob
    
    message(cat("Starting theta posterior probability: ",posterior,sep="\t"))
    
####################
    ## PRE ALLOCATE MEMORY
####################
    ## Create empty chain to store every iteration for the adaptive period and burnin
    #opt_chain <- matrix(nrow=burnin + adaptive_period,ncol=unfixed_par_length)
    #opt_chain <- rep(list(opt_chain),length(temperatures))
    chain_index <- 1
    
    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=param_length+4)

## Set up initial csv file
    ## Store posterior (called lnlike), likelihood ad prior
    chain_colnames <- c("sampno",par_names,"lnlike","likelihood","prior_prob")
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,current_pars,posterior,total_likelihood, prior_prob)
    colnames(tmp_table) <- chain_colnames
    
    ## Write starting conditions to file
    data.table::fwrite(as.data.frame(tmp_table),file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)
    save_infection_history_to_disk(infection_histories, infectionHistory_file, 1, append=FALSE,colNames=TRUE)

    ## Initial indexing parameters
    no_recorded <- 1
    no_recorded_infection_history <- 1
    sampno <- 2
    par_i <- 1
    i <- 1
    chain_index <- 1
    last_index <- 1
    offset <- 0
    potential_swaps <- swaps <- 0

    ## set seed
    if(!missing(seed)){
        ## using an integer
        if(length(seed) == 1){
            set.seed(seed)
            ## or a previous state
        } else {
            .Random.seed <- seed
        }
    }
  

    
    
    ## initialise MCMC
    mcmc_list <- list("i"=i,"par_i" = par_i, "current_pars" = current_pars,
                      "infection_histories"=infection_histories,
                      "likelihoods"=likelihoods,"total_likelihood" = total_likelihood,
                      "prior_prob"=prior_prob,"posterior"=posterior,
                      "tempaccepted" = tempaccepted,"tempiter" = tempiter,
                      "steps"=NULL,"temp"=1)
    ## replicate list for parallel tempering
    mcmc_list <- rep(list(mcmc_list),length(temperatures))

    ## start values for parallel tempering
    mcmc_list <- Map(function(x,y) modifyList(x,list(current_pars = y)), mcmc_list, start_pars)
    mcmc_list <- Map(function(x,y) modifyList(x,list(steps = y)), mcmc_list, steps)
    mcmc_list <- Map(function(x,y) modifyList(x,list(temp = y)), mcmc_list, temperatures)
    ##mcmc_list <- Map(function(x,y) modifyList(x,list(infection_histories = y)), mcmc_list, infection_histories)
    
    
    

    
                                        # registerDoMC(cores=6)
    
    run_multiple_iter <- function(tmp_list){
        startI <- tmp_list$i
        for(i in startI:(startI+n_iter)){
            tmp_list <- run_MCMC_single_iter(i, tmp_list$par_i,
                                             tmp_list$current_pars, tmp_list$infection_histories,
                                             tmp_list$likelihoods, tmp_list$total_likelihood,
                                             tmp_list$prior_prob, tmp_list$posterior,
                                             tmp_list$tempaccepted, tmp_list$tempiter,
                                             tmp_list$steps, tmp_list$temp)
            startI <- startI + 1
        }
        tmp_list$i <- startI
        tmp_list        
    }
    
    run_MCMC_single_iter <- lapply(1:length(temperatures),
                                   function(x) create_run_MCMC_single_iter_fn_fast(unfixed_pars,unfixed_par_length,
                                                                                   lower_bounds,upper_bounds,
                                                                                   age_mask, strain_mask,
                                                                                   posterior_simp,
                                                                                   extra_probabilities,
                                                                                   proposal_gibbs,
                                                                                   alpha,beta,histSampleProb,
                                                                                   nInfs_vec,
                                                                                   swapPropn,moveSize,n_alive,
                                                                                   switch_sample,
                                                                                   hist_switch_prob,
                                                                                   year_swap_propn))
    #cl <- makeForkCluster(6, useXDR=FALSE)

    ## main body of running MCMC
    while (i <= (iterations+adaptive_period)){
        for(jh in 1:length(mcmc_list)) mcmc_list[[jh]][["i"]] <- i
##############
        ## Can probably parallelise this bit
        mcmc_list <- Map(do.call, run_MCMC_single_iter, mcmc_list)
                                        #mcmc_list <- lapply(mcmc_list, run_multiple_iter)
                                        #mcmc_list <- lapply(mcmc_list, do.call, run_MCMC_single_iter[[1]])
                                        #mcmc_list <- mclapply(mcmc_list, run_multiple_iter)
                                        #mcmc_list <- parLapply(cl, mcmc_list, run_multiple_iter)
                                        #foreach(jh=1:length(mcmc_list)) %dopar%{
                                        #    mcmc_list[[jh]] <- run_multiple_iter(mcmc_list[[jh]])
                                        #}
                                        #mcmc_list <- clusterApply(NULL, mcmc_list, run_multiple_iter)
                                        #mcmc_list <- Map(do.call, run_MCMC_single_iter, mcmc_list)
##############
        ## perform parallel tempering
        if(i %% parallel_tempering_iter == 0){
            parallel_tempering_list <- parallel_tempering(mcmc_list, temperatures, offset)
            mcmc_list <- parallel_tempering_list$mcmc_list
            swaps <- swaps + parallel_tempering_list$swaps
            potential_swaps <- potential_swaps + 1
            offset <- 1 - offset
        }
##############################
        ## SAVE STEP
##############################
        ## If current iteration matches with recording frequency, store in the chain.
        ## If we are at the limit of the save block,
        ## save this block of chain to file and reset chain
        ## Save theta
        if(i %% thin ==0){
            current_pars <- mcmc_list[[1]][["current_pars"]]
            total_likelihood <- mcmc_list[[1]][["total_likelihood"]]
            prior_prob <- mcmc_list[[1]][["prior_prob"]]
            posterior <- mcmc_list[[1]][["posterior"]]
            save_chain[no_recorded,1] <- sampno
            save_chain[no_recorded,2:(ncol(save_chain)-3)] <- current_pars
            save_chain[no_recorded,ncol(save_chain)-2] <- posterior
            save_chain[no_recorded,ncol(save_chain)-1] <- total_likelihood
            save_chain[no_recorded,ncol(save_chain)] <- prior_prob
            no_recorded <- no_recorded + 1
        }

        ## If within adaptive period, need to do some adapting!
        if(i <= adaptive_period){
            ## Save each step
            #for(jh in 1:length(opt_chain)) opt_chain[[jh]][chain_index,] <- mcmc_list[[jh]][["current_pars"]][unfixed_pars]
            ## If in an adaptive step
            if(chain_index %% opt_freq == 0){
                reset_acceptance <- function(mcmc_list, reset){
                    mcmc_list[["tempaccepted"]] <- mcmc_list[["tempiter"]] <- reset
                    mcmc_list
                }
                ## Current acceptance rate
                pcur <- lapply(mcmc_list, function(x) x[["tempaccepted"]] / x[["tempiter"]])
                scale_univariate <- function(steps, popt, pcur, unfixed_pars){
                    steps[unfixed_pars] <- vapply(unfixed_pars,function(x) scaletuning(steps[x],popt,pcur[x]),
                                                  double(1))
                    steps
                }
                steps <- Map(function(x,y) scale_univariate(x, popt, y, unfixed_pars), steps, pcur)

                mcmc_list <- Map(function(x,y) modifyList(x,list(steps = y)), mcmc_list, steps)
                
                message(cat("Pcur: ", pcur[[1]][unfixed_pars],sep="\t"))
                message(cat("Step sizes: ", steps[[1]][unfixed_pars],sep="\t"))
                
                mcmc_list <- lapply(mcmc_list, function(x) reset_acceptance(x, reset))
                ## calibrate temperatures
                swap_ratio <- swaps / potential_swaps
                message(cat("Swap ratio: ", swap_ratio,sep="\t"))
                temperatures <- calibrate_temperatures(temperatures, swap_ratio)
                message(cat("Temperatures: ", temperatures,sep="\t"))
                mcmc_list <- Map(function(x,y) modifyList(x,list(temp = y)), mcmc_list, temperatures)
                swaps <- potential_swaps <- 0
            }
            chain_index <- chain_index + 1
        }
        if(i %% save_block == 0) message(cat("Current iteration: ", i, sep="\t"))
        if(no_recorded == save_block){
            data.table::fwrite(as.data.frame(save_chain[1:(no_recorded-1),]),file=mcmc_chain_file,
                               col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            save_chain <- empty_save_chain
            no_recorded <- 1
        }
        
        ## Save infection histories
        if(i %% histTabThin == 0){
            infection_histories <- mcmc_list[[1]][["infection_histories"]]
            save_infection_history_to_disk(infection_histories, infectionHistory_file, sampno)                                      
        }
        sampno <- sampno + n_iter
        i <- i + 1

    }
    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        data.table::fwrite(as.data.frame(save_chain[1:(no_recorded-1),]),file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }
    p_accept <- mcmc_list[[1]][["tempaccepted"]] / mcmc_list[[1]][["tempiter"]]
    p_accept <- p_accept[unfixed_pars]
    
    current_pars <- lapply(mcmc_list, function(x)x$current_pars)
    
    #parallelStop()
    
    return(list("chain_file" = mcmc_chain_file,
                "history_file"=infectionHistory_file,
                "current_pars" = current_pars,
                "steps" = steps, 
                "p_accept" = p_accept,
                "temperatures" = temperatures,
                "seed" = .Random.seed))
}



