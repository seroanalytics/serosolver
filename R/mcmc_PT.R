run_MCMC_pt <- function(parTab,
                        data,
                        antigenicMap,
                        mcmcPars=c(),
                        startInfHist=NULL,
                        filename="test",
                        CREATE_POSTERIOR_FUNC,
                        PRIOR_FUNC=NULL,
                        version=2,
                        seed,
                        mu_indices=NULL,
                        measurement_indices=NULL,
                        measurement_random_effects=FALSE,
                        OPT_TUNING=0.1,
                        temp=1,
                        solve_likelihood=TRUE,
                        n_alive=NULL,
                        ...){
    ## Extract MCMC parameters
    mcmcPars_used <- c("iterations"=50000,"popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,"thin"=1,
                       "adaptive_period"=10000,
                       "save_block"=100,"thin2"=10,"histSampleProb"=1,"switch_sample"=2, "burnin"=0, 
                       "inf_propn"=1, "moveSize"=5,"histOpt"=1,"swapPropn"=0.5,
                       "hist_switch_prob"=0,"year_swap_propn"=1)
    mcmcPars_used[names(mcmcPars)] <- mcmcPars

    ## Extract MCMC parameters
    iterations <- mcmcPars_used["iterations"] # How many iterations to run after adaptive period
    popt <- mcmcPars_used["popt"] # Desired optimal acceptance rate
    popt_hist <- mcmcPars_used["popt_hist"]
    opt_freq<- mcmcPars_used["opt_freq"] # How often to adjust step size
    thin <- mcmcPars_used["thin"] # Save only every nth iterations for theta sampling
    adaptive_period<- mcmcPars_used["adaptive_period"] # How many iterations for adaptive period
    save_block <- mcmcPars_used["save_block"] # How many post-thinning iterations to store before saving to disk
    histTabThin <- mcmcPars_used["thin2"] # Save only every nth iterations for infection history sampling
    histSampleProb <- mcmcPars_used["histSampleProb"] # What proportion of infection histories to sample each step
    switch_sample <- mcmcPars_used["switch_sample"] # Resample infection histories every n iterations
    burnin <- mcmcPars_used["burnin"] # Run this many iterations before attempting adaptation. Idea is to reduce getting stuck in local maxima
    moveSize <- mcmcPars_used["moveSize"] # Number of infections to move/remove/add in each proposal step
    inf_propn <- mcmcPars_used["inf_propn"] # Number of infections to move/remove/add in each proposal step
    nInfs <- floor(length(antigenicMap$inf_years)*inf_propn)
    histOpt <- mcmcPars_used["histOpt"] # Should infection history proposal step be adaptive?
    swapPropn <- mcmcPars_used["swapPropn"] # If using gibbs, what proportion of proposals should be swap steps?
    hist_switch_prob <- mcmcPars_used["hist_switch_prob"] # If using gibbs, what proportion of iterations should be swapping contents of two time periods?
    year_swap_propn <- mcmcPars_used["year_swap_propn"] # If gibbs and swapping contents, what proportion of these time periods should be swapped?
    temperatures <- mcmcPars["temperature"]
    parallel_tempering_iter <- mcmcPars["parallel_tempering_iter"]
###################################################################

    current_pars <- parTab[[1]]$values
    start_pars <- lapply(parTab, function(x) x$values)
    steps <- lapply(parTab, function(x) x$steps)
    parTab <- parTab[[1]]

    param_length <- nrow(parTab)
    unfixed_pars <- which(parTab$fixed == 0)
    unfixed_par_length <- nrow(parTab[parTab$fixed== 0,])
    par_names <- as.character(parTab$names)
    
    ## Parameter constraints
    lower_bounds <- parTab$lower_bound
    upper_bounds <- parTab$upper_bound
    fixed <- parTab$fixed

    ## If using lambda terms, pull their indices out of the parameter table
    lambda_indices <- NULL
    if("lambda" %in% par_names){
        lambda_indices <- which(parTab$names == "lambda")
    }
    alpha <- parTab[parTab$names == "alpha","values"]
    beta <- parTab[parTab$names == "beta","values"]
    
    tempaccepted <- tempiter <- integer(param_length)
    reset <- integer(param_length)
    reset[] <- 0

    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste0(filename,"_chain.csv")
    infectionHistory_file <- paste0(filename,"_infectionHistories.csv")


    strainIsolationTimes <- unique(antigenicMap$inf_years) # How many strains are we testing against and what time did they circulate
    samplingTimes <- unique(data$sample) # What are the range of sampling times?
    n_strain <- length(strainIsolationTimes) # How many strains could an individual see?
    n_indiv <- length(unique(data$individual)) # How many individuals in the data?
    n_groups <- length(unique(data$group)) # How many groups in the data?
    individuals <- 1:n_indiv # Create vector of individuals
    groups <- 1:n_groups # Create vector of groups

###############
    ## Create age mask
    ## -----------------------
    ## Note that ages for all groups must be from same reference point
    ## -----------------------
###############
    if(!is.null(titreDat$DOB)){
        DOBs <- unique(titreDat[,c("individual","DOB")])[,2]
    } else {
        DOBs <- rep(min(strainIsolationTimes), n_indiv)
    }
    ageMask <- create_age_mask(DOBs, strainIsolationTimes)
    ## Create strain mask
    strainMask <- create_strain_mask(titreDat,strainIsolationTimes)
    masks <- data.frame(cbind(ageMask, strainMask))
    ## Number of people that were born before each year and have had a sample taken since that year happened
    if(is.null(n_alive)) n_alive <- sapply(seq(1,length(strainIsolationTimes)), function(x) nrow(masks[masks$ageMask <=x & masks$strainMask >= x,]))   

    
    ## Create posterior calculating function
    posterior_simp <- protect(CREATE_POSTERIOR_FUNC(parTab,
                                                    data,
                                                    antigenicMap,
                                                    version,
                                                    solve_likelihood,
                                                    ageMask,
                                                    measurement_indices=measurement_indices,
                                                    mu_indices=mu_indices,
                                                    n_alive=n_alive,
                                                    function_type=1,
                                                    ...))
    if(!is.null(CREATE_PRIOR_FUNC)){
        prior_func <- CREATE_PRIOR_FUNC(parTab)
    }
    ## If using gibbs proposal on infHist, create here
    proposal_gibbs <- protect(CREATE_POSTERIOR_FUNC(parTab,
                                                        titreDat,
                                                        antigenicMap,
                                                        version,
                                                        solve_likelihood,
                                                        ageMask,
                                                        measurement_indices_by_time=measurement_indices,
                                                        mu_indices=mu_indices,
                                                        n_alive=n_alive,
                                                        function_type=2,
                                                        ...))
    
   ## If using random effects on mu, need to include hyperprior term on mu
    ## We can't do this in the main posterior function, because this term
    ## applies to the overall posterior whereas the main posterior function
    ## returns each individual's posterior
    if (!is.null(mu_indices)) {
        prior_mu <- create_prior_mu(parTab)
    }
    if (measurement_random_effects) {
        prior_shifts <- create_prob_shifts(parTab)
    }
######################
    
    ## Setup initial conditions
    infectionHistories = startInfHist
    if (is.null(startInfHist)) {
        infectionHistories <- setup_infection_histories_new(titreDat, strainIsolationTimes, space=5,titre_cutoff=3)
    }

    ## Initial likelihood
    likelihoods <- posterior_simp(current_pars,infectionHistories)/temp
    ## Initial total likelihood
    total_likelihood <- sum(likelihoods)
    n_alive_tot <- sum(n_alive)
    ## Create closure to add extra prior probabilities, to avoid re-typing later
    extra_probabilities <- function(prior_pars, prior_infHist){
        prior_probab <- 0
        if(histProposal == 2){
            if(prior_on_total){
                prior_probab <- prior_probab + inf_mat_prior_total_cpp(prior_infHist, n_alive_tot, alpha, beta)
            } else {
                prior_probab <- prior_probab + inf_mat_prior_cpp(prior_infHist, n_alive, alpha, beta)
            }
        }
        if(!is.null(CREATE_PRIOR_FUNC)) prior_probab <- prior_probab + prior_func(prior_pars)
        if(!is.null(mu_indices)) prior_probab <- prior_probab + prior_mu(prior_pars)
        if(measurement_random_effects) prior_probab <- prior_probab + prior_shifts(prior_pars)
        prior_probab
    }
    ## Initial total prior prob
    prior_prob <- extra_probabilities(current_pars, infectionHistories)

    ## Initial posterior prob
    posterior <- total_likelihood + prior_prob
    
    message(cat("Starting theta posterior probability: ",posterior,sep="\t"))
    
####################
    ## PRE ALLOCATE MEMORY
####################
    ## Create empty chain to store every iteration for the adaptive period and burnin
    opt_chain <- matrix(nrow=burnin + adaptive_period,ncol=unfixed_par_length)
    opt_chain <- rep(list(opt_chain),length(temperatures))
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
    save_infHist_to_disk(infectionHistories, infectionHistory_file, 1, append=FALSE,colNames=TRUE)

    ## Initial indexing parameters
    no_recorded <- 1
    no_recorded_infHist <- 1
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

    run_MCMC_single_iter <- lapply(seq_along(temperatures),
                                   function(x) create_run_MCMC_single_iter_fn(unfixed_pars,unfixed_par_length,
                                                                              lower_bounds,upper_bounds,
                                                                              ageMask, strainMask,
                                                                              posterior_simp, extra_probabilities,
                                                                              proposal_gibbs,
                                                                              alpha,beta,histSampleProb,nInfs_vec,
                                                                              swapPropn,moveSize,n_alive,
                                                                              switch_sample,
                                                                              hist_switch_prob,
                                                                              year_swap_propn))
    ## initialise MCMC
    mcmc_list <- list("i"=i,"par_i" = par_i, "current_pars" = current_pars,
                      "infectionHistories"=infectionHistories,
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

    ## main body of running MCMC
    while (i <= (iterations+adaptive_period)){
        for(jh in 1:length(mcmc_list)) mcmc_list[[jh]][["i"]] <- i
        mcmc_list <- Map(do.call, run_MCMC_single_iter, mcmc_list)

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
            likelihood <- mcmc_list[[1]][["likelihood"]]
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
            for(jh in 1:length(opt_chain)) opt_chain[[jh]][chain_index,] <- mcmc_list[[jh]][["current_pars"]][unfixed_pars]
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
            infectionHistories <- mcmc_list[[1]][["infectionHistories"]]
            save_infHist_to_disk(infectionHistories, infectionHistory_file, sampno)                                      
        }
        sampno <- sampno + 1
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
    return(list("chain_file" = mcmc_chain_file,
                "history_file"=infectionHistory_file,
                "current_pars" = current_pars,
                "steps" = steps, 
                "p_accept" = p_accept,
                "temperatures" = temperatures,
                "seed" = .Random.seed))
}



