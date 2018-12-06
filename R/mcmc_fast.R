#' Adaptive Metropolis-within-Gibbs/Metropolis Hastings Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. The algorithm then uses \code{\link{univ_proposal}} or \code{\link{mvr_proposal}} to explore parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename. This version of the algorithm is also designed to explore posterior densities for infection histories
#' @param par_tab The parameter table controlling information such as bounds, initial values etc
#' @param data The data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time); run (integer giving the repeated number of this titre); DOB (integer giving date of birth matching time units used in model)
#' @param antigenic_map A data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years 
#' @param mcmcPars Named vector named vector with parameters for the MCMC procedure. See details
#' @param startInfHist Infection history matrix to start MCMC at. Can be left NULL
#' @param filename The full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param CREATE_POSTERIOR_FUNC Pointer to posterior function used to calculate a likelihood. This will probably be \code{\link{create_posterior_func}}, but if using random effects on mu will use \code{\link{create_post_func_mu}}
#' @param PRIOR_FUNC User function of prior for model parameters. Should take parameter values only
#' @param version Which version of the posterior function to use? See \code{\link{create_post_func}}
#' @param mu_indices For random effects on boosting parameter, mu. Vector of indices of length equal to number of circulation times. If random mus are included in the parameter table, this vector specifies which mu to use for each circulation year. For example, if years 1970-1976 have unique boosting, then mu_indices should be c(1,2,3,4,5,6). If every 3 year block shares has a unique boosting parameter, then this should be c(1,1,1,2,2,2)
#' @param measurement_indices For measurement bias function. Vector of indices of length equal to number of circulation times. For each year, gives the index of parameters named "rho" that correspond to each time period
#' @param measurement_random_effects Boolean indicating if measurement bias is a random effects term. If TRUE adds a component to the posterior calculation that calculates the probability of a set of measurement shifts "rho", given a mean and standard deviation
#' @param OPT_TUNING Constant describing the amount of leeway when adapting the proposals steps to reach a desired acceptance rate (ie. does not change step size if within OPT_TUNING of the specified acceptance rate)
#' @param temp Temperature term for parallel tempering, raises likelihood to this value. Just used for testing at this point
#' @param solve_likelihood if FALSE, returns only the prior and does not solve the likelihood
#' @param ... Other arguments to pass to CREATE_POSTERIOR_FUNC
#' @return A list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) a full file path at which the infection history chain is saved as a .csv file; 3) the last used covariance matrix; 4) the last used scale/step size (if multivariate proposals)
#' @details
#' The `mcmcPars` argument has the following options:
#'  * iterations (number of post adaptive period iterations to run)
#'  * adaptive_period (for this many iterations, change proposal step size adaptively every `opt_freq` iterations)
#'  * opt_freq (adapt proposal step size every opt_freq iterations)
#'  * thin (save every n iterations)
#'  * thin_hist (infection history thinning)  
#'  * save_block (after this many iterations (post thinning), save to disk)
#'  * popt (desired acceptance rate for parameters)
#'  * popt_hist (desired acceptance rate for infection histories)
#'  * switch_sample (resample inf histories every n iterations)
#'  * inf_propn (proportion of infection times to resample for each individual at each iteration)
#'  * moveSizes (number of infection years/months to move when performing infection history swap step)
#'  * histProposal (which infection history proposal version to use, see \code{\link{describe_proposals}}
#'  * histOpt (if 1, performs adaptive infection history proposals. If 0, retains the starting infection history proposal parameters
#'  * swapPropn (if using gibbs sampling of infection histories, what proportion of proposals should be swap steps)
#'  * hist_switch_prob (proportion of infection history proposal steps to swap year_swap_propn of two time periods' contents)
#'  * year_swap_propn (when swapping contents of two time points, what proportion of individuals should have their contents swapped)
#' @md
#' @export
run_MCMC_fast <- function(par_tab,
                          titre_dat,
                          antigenic_map,
                          mcmcPars=c(),
                          startInfHist=NULL,
                          filename="test",
                          CREATE_POSTERIOR_FUNC=create_posterior_func_fast,
                          CREATE_PRIOR_FUNC=NULL,
                          mu_indices=NULL,
                          measurement_indices=NULL,
                          measurement_random_effects=FALSE,
                          OPT_TUNING=0.1,
                          temp=1,
                          solve_likelihood=TRUE,
                          n_alive=NULL,
                          message_slack=FALSE,
                          message_slack_pars=NULL,
                          ...){
    ## If we want to message ourselves on slack with intermediate progress
    if(message_slack && !is.null(message_slack_pars)){
        slackr::slackrSetup(api_token=message_slack_pars$api_token,
                            channel=message_slack_pars$channel,
                            username=message_slack_pars$username)        
    }

    
    
    ## Sort out MCMC parameters --------------------------------------
    ###################################################################
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
    nInfs <- floor(length(antigenic_map$inf_years)*inf_propn)
    histOpt <- mcmcPars_used["histOpt"] # Should infection history proposal step be adaptive?
    swapPropn <- mcmcPars_used["swapPropn"] # If using gibbs, what proportion of proposals should be swap steps?
    hist_switch_prob <- mcmcPars_used["hist_switch_prob"] # If using gibbs, what proportion of iterations should be swapping contents of two time periods?
    year_swap_propn <- mcmcPars_used["year_swap_propn"] # If gibbs and swapping contents, what proportion of these time periods should be swapped?
    ###################################################################    

    ## Sort out which version to run --------------------------------------
    propPrint <- "Using integrated FOI prior on infection history, with gibbs sampling of infections"
    histProposal <- 2
    
    message(propPrint)
    
    ## Extract parameter settings
    par_names <- as.character(par_tab$names) # Parameter names
    
    param_length <- nrow(par_tab)
    unfixed_pars <- which(par_tab$fixed == 0) # Indices of free parameters
    unfixed_par_length <- nrow(par_tab[par_tab$fixed== 0,]) # How many free parameters?
    current_pars <- par_tab$values # Starting parameters
    ## Parameter constraints
    lower_bounds <- par_tab$lower_bound # Parameters cannot step below this
    upper_bounds <- par_tab$upper_bound # Parameters cannot step above this
    steps <- par_tab$steps # How far to step on unit scale to begin with?
    fixed <- par_tab$fixed # Index which parameters are fixed
    
    alpha <- par_tab[par_tab$names == "alpha","values"]
    beta <- par_tab[par_tab$names == "beta","values"]

    ## To store acceptance rate of entire time period infection history swaps
    infection_historySwapN <- infection_historySwapAccept <- 0
    ## Arrays to store acceptance rates
    ## If univariate proposals, store vector of acceptances
    tempaccepted <- tempiter <- integer(param_length)
    reset <- integer(param_length)
    reset[] <- 0
    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste0(filename,"_chain.csv")
    infectionHistory_file <- paste0(filename,"_infection_histories.csv")
    
###############
    ## Extract titre_dat parameters
##############
    ## Check the titre_dat input 
    check_data(titre_dat)
    
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
    histiter <- histaccepted <- histiter_add <- histaccepted_add <- histiter_move <- histaccepted_move <- histreset <- integer(n_indiv)
    
    nInfs_vec <- rep(nInfs, n_indiv) # How many infection history moves to make with each proposal
    moveSizes <- rep(moveSize, n_indiv) # How many years to move in smart proposal step
###############
    ## Create age mask
    ## -----------------------
    ## Note that DOBs for all groups must be from same reference point
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
    proposal_gibbs <- protect(CREATE_POSTERIOR_FUNC(par_tab,
                                                    titre_dat,
                                                    antigenic_map,
                                                    version,
                                                    solve_likelihood,
                                                    age_mask,
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
        prior_mu <- create_prior_mu(par_tab)
    }
    if (measurement_random_effects) {
        prior_shifts <- create_prob_shifts(par_tab)
    }
######################
    
    ## Setup initial conditions
    infection_histories = startInfHist
    if (is.null(startInfHist)) {
        infection_histories <- setup_infection_histories_new(titre_dat, strain_isolation_times, space=5,titre_cutoff=3)
    }
   
    ## Initial likelihood
    likelihoods <- posterior_simp(current_pars,infection_histories)/temp
    ## Initial total likelihood
    total_likelihood <- sum(likelihoods)
    n_alive_tot <- sum(n_alive)
    ## Create closure to add extra prior probabilities, to avoid re-typing later
    extra_probabilities <- function(prior_pars, prior_infection_history){
        prior_probab <- 0
        prior_probab <- prior_probab + inf_mat_prior_cpp(prior_infection_history, n_alive, alpha, beta)
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
###############
    
####################
    ## PRE ALLOCATE MEMORY
####################
    ## Create empty chain to store every iteration for the adaptive period and burnin
    opt_chain <- matrix(nrow=burnin + adaptive_period,ncol=unfixed_par_length)

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
    chain_index <- 1
    last_index <- 1

    cov_mat0 <- diag(unfixed_pars)
#####################
    ## MCMC ALGORITHM
#####################
    log_prob <- 0
    for(i in 1:(iterations + adaptive_period + burnin)){
        ## Whether to swap entire year contents or not - only applies to gibbs sampling
        histSwitchProb <- runif(1)
        if(i %% save_block == 0) message(cat("Current iteration: ", i, sep="\t"))
        if(message_slack && i %% message_slack_pars$message_freq == 0){
            text_slackr(paste0(message_slack_pars$username, " iteration ", i),channel=message_slack_pars$channel,
                        username=message_slack_pars$username)
        }
######################
        ## PROPOSALS
######################
        ## If updating theta
        if (i %% switch_sample == 0) {
            ## If all pars are fixed 
            if (length(unfixed_pars)==0) {
                proposal <- current_pars ## Set proposal to be current parameters
                tempiter <- tempiter + 1
            } else { ##Else propose parameters
                ## For each parameter (Gibbs)
                j <- unfixed_pars[par_i]
                par_i <- par_i + 1
                if(par_i > unfixed_par_length) par_i <- 1
                proposal <- univ_proposal(current_pars, lower_bounds, upper_bounds, steps,j) 
                tempiter[j] <- tempiter[j] + 1
                ## If using multivariate proposals
            }
            ## Calculate new likelihood for these parameters
            new_likelihoods <- posterior_simp(proposal,infection_histories)/temp # For each individual
            new_total_likelihood <- sum(new_likelihoods) # Total
            new_prior_prob <- extra_probabilities(proposal, infection_histories) # Prior
            new_posterior <- new_total_likelihood + new_prior_prob # Posterior
            ## Otherwise, resample infection history
        } else {
            ## Choose a random subset of individuals to update
            indivSubSample <- sample(1:n_indiv, ceiling(histSampleProb*n_indiv))
            ## Need to temporarily store current parameters as new pars, as
            ## might change with lambda swap step
            proposal <- current_pars            
            
            if(histSwitchProb > hist_switch_prob){
                prop_gibbs <- proposal_gibbs(proposal, infection_histories,
                                             likelihoods,
                                             indivSubSample, 
                                             alpha, beta,
                                             nInfs_vec,swapPropn,moveSize,
                                             temp)
                new_likelihoods <- prop_gibbs$old_probs
                newInfectionHistories <- prop_gibbs$new_infection_history
            } else {
                newInfectionHistories <- inf_hist_swap(infection_histories, age_mask, strain_mask,
                                                       year_swap_propn, moveSize)
                new_likelihoods <- posterior_simp(proposal, newInfectionHistories)/temp
                if(!identical(newInfectionHistories, infection_histories)) infection_historySwapN <- infection_historySwapN + 1
            }
            
            ## Calculate new likelihood with these infection histories
            new_total_likelihood <- sum(new_likelihoods)
            new_prior_prob <- extra_probabilities(proposal, newInfectionHistories)
            new_posterior <- new_total_likelihood + new_prior_prob
        }

#############################
        ## METROPOLIS HASTINGS STEP
        #############################
        ## Check that all proposed parameters are in allowable range
        ## Skip if any parameters are outside of the allowable range
        if(i %% switch_sample == 0){
            log_prob <- new_posterior-posterior
            if(!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)){
                log_prob <- min(log_prob, 0)
                if(log(runif(1)) < log_prob){
                    if(!any(proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                            proposal[unfixed_pars] > upper_bounds[unfixed_pars])){
                        ## Accept with probability 1 if better, or proportional to
                        ## difference if not
                        current_pars <- proposal
                        ## Store acceptances
                        ## If all parameters are fixed, then we 'accept'
                        if(length(unfixed_pars)==0){ 
                            tempaccepted <- tempaccepted + 1
                        }else{
                            tempaccepted[j] <- tempaccepted[j] + 1
                        }
                        likelihoods <- new_likelihoods
                        prior_prob <- new_prior_prob
                        total_likelihood <- new_total_likelihood
                        posterior <- new_posterior
                    }
                }
            }
        } else {
            if(histSwitchProb > hist_switch_prob){
                ## MH step for each individual
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
                            if(!any(proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                                    proposal[unfixed_pars] > upper_bounds[unfixed_pars])){
                                infection_historySwapAccept <- infection_historySwapAccept + 1
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
            }
        }
        
        
##############################
        ## SAVE STEP
##############################
        ## If current iteration matches with recording frequency, store in the chain.
        ## If we are at the limit of the save block,
        ## save this block of chain to file and reset chain
        
        ## Save theta
        if(i %% thin ==0){
            save_chain[no_recorded,1] <- sampno
            save_chain[no_recorded,2:(ncol(save_chain)-3)] <- current_pars
            save_chain[no_recorded,ncol(save_chain)-2] <- posterior
            save_chain[no_recorded,ncol(save_chain)-1] <- total_likelihood
            save_chain[no_recorded,ncol(save_chain)] <- prior_prob
            no_recorded <- no_recorded + 1
        }

        ## Save infection histories
        if(i %% histTabThin == 0){
            save_infection_history_to_disk(infection_histories, infectionHistory_file, sampno)
            ##historyTab[no_recorded_infection_history:(no_recorded_infection_history + n_indiv-1),1:n_strain] <- infection_histories
            ##historyTab[no_recorded_infection_history:(no_recorded_infection_history + n_indiv-1),n_strain+1] <- individuals
            ##historyTab[no_recorded_infection_history:(no_recorded_infection_history + n_indiv-1),n_strain+2] <- sampno
            ##no_recorded_infection_history <- no_recorded_infection_history + n_indiv
        }
        
##############################
        ## ADAPTIVE PERIOD
##############################
        ## If within adaptive period, need to do some adapting!
        if(i > (adaptive_period + burnin) & i %% opt_freq == 0){
            pcur <- tempaccepted/tempiter ## get current acceptance rate
            message(cat("Pcur: ", signif(pcur,3),sep="\t"))
            message(cat("Step sizes: ", signif(steps,3),sep="\t"))
            message(cat("Inf hist swap pcur: ", signif(infection_historySwapAccept/infection_historySwapN,3),sep="\t"))
            infection_historySwapAccept <- infection_historySwapN <- 0
            tempaccepted <- tempiter <- reset

            ## Have a look at the acceptance rates for infection histories
            pcurHist <- histaccepted/histiter ## Overall
            pcurHist_add <- histaccepted_add/histiter_add ## For adding
            pcurHist_move <- histaccepted_move/histiter_move ## For moving
            histiter <- histaccepted <- histaccepted_add <- histaccepted_move <-
                histiter_add <- histiter_move <- histreset
        }
        if(i > burnin & i <= (adaptive_period + burnin)){
            ## Current acceptance rate
            pcur <- tempaccepted/tempiter
            ## Save each step
            opt_chain[chain_index,] <- current_pars[unfixed_pars]
            ## If in an adaptive step
            if(chain_index %% opt_freq == 0){
                ## If using univariate proposals
                ## For each non fixed parameter, scale the step size
                for(x in unfixed_pars) steps[x] <- scaletuning(steps[x],popt,pcur[x])

                pcurHist <- histaccepted/histiter
                pcurHist_add <- histaccepted_add/histiter_add
                pcurHist_move <- histaccepted_move/histiter_move

                ## NOTE THAT THIS IS ONLY RELEVANT TO INFECTION HISTORY PROPOSAL 1 & 3
                if(histProposal != 2){
                    #message(cat("Hist acceptance add: ", signif(pcurHist_add,3), cat="\t"))
                    #message(cat("Hist acceptance move: ", signif(pcurHist_move,3), cat="\t"))
                    message(cat("Mean hist acceptance: ", signif(mean(pcurHist),3),cat="\t"))
                    histiter <- histaccepted <- histaccepted_add <- histaccepted_move <- histiter_add <- histiter_move <- histreset
                }

                ## If adaptive infection history proposal
                if(histOpt == 1){
                    ## Increase or decrease the number of infection history locations
                    ## being changed to modify acceptance rate. If not accepting enough,
                    ## reduce number. If accepting too many, increase number                    
                    nInfs_vec[which(pcurHist_add < popt_hist*(1-OPT_TUNING))] <- nInfs_vec[which(pcurHist_add< popt_hist*(1-OPT_TUNING))] - 1
                    nInfs_vec[which(pcurHist >= popt_hist*(1+OPT_TUNING))] <- nInfs_vec[which(pcurHist >= popt_hist*(1+OPT_TUNING))] +1
                    nInfs_vec[nInfs_vec < 1] <- 1
                    #moveSizes[which(pcurHist_move< popt_hist*(1-OPT_TUNING))] <- moveSizes[which(pcurHist_move < popt_hist*(1-OPT_TUNING))] - 1
                    #moveSizes[which(pcurHist >= popt_hist*(1+OPT_TUNING))] <- moveSizes[which(pcurHist >= popt_hist*(1+OPT_TUNING))] +1
                    #moveSizes[moveSizes < 1] <- 1

                    for(ii in seq_along(nInfs_vec)){
                        moveSizes[ii] <- min(moveSizes[ii], length(age_mask[ii]:strain_mask[ii]))
                        nInfs_vec[ii] <- min(nInfs_vec[ii],length(age_mask[ii]:strain_mask[ii]))
                    }
                }
                message(cat("nInfs: ", head(nInfs_vec), sep="\t"))
                message(cat("Move sizes: ", head(moveSizes), sep="\t"))
                message(cat("Inf hist swap pcur: ", signif(infection_historySwapAccept/infection_historySwapN,3),sep="\t"))
                message(cat("Pcur: ", signif(pcur,3),sep="\t"))
                message(cat("Step sizes: ", signif(steps,3),sep="\t"))
                if(all(pcur == 0)){
                    message("Warning: acceptance rates are 0. Might be an error with the theta proposal?")
                    if(message_slack){
                        text_slackr(paste0("Warning in ", message_slack_pars$username,": acceptance rates are 0. Might be an error with the theta proposal?"),
                                    channel=message_slack_pars$channel,
                                    username=message_slack_pars$username)
                    }
                }
                
                tempaccepted <- tempiter <- reset
            }
            chain_index <- chain_index + 1
        }
#######################
        ## HOUSEKEEPING
#######################
        if(no_recorded == save_block){
            data.table::fwrite(as.data.frame(save_chain[1:(no_recorded-1),]),file=mcmc_chain_file,
                               col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            save_chain <- empty_save_chain
            no_recorded <- 1
        }
                                        #if((no_recorded_infection_history-1)/n_indiv == save_block){
                                        #    data.table::fwrite(as.data.frame(historyTab[1:(no_recorded_infection_history-1),]), file=infectionHistory_file,
                                        #                col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
                                        #    historyTab <- emptyHistoryTab
                                        #    no_recorded_infection_history <- 1
                                        #}
        sampno <- sampno + 1
    }

    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        data.table::fwrite(as.data.frame(save_chain[1:(no_recorded-1),]),file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }
                                        #if(no_recorded_infection_history > 2){
                                        #    data.table::fwrite(as.data.frame(historyTab), file=infectionHistory_file,
                                        #                col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
                                        #    historyTab <- emptyHistoryTab
                                        #    no_recorded_infection_history <- 1
                                        #}

    if(message_slack && !is.null(message_slack_pars)){
        text_slackr(paste0(message_slack_pars$username, " is done!"))
    }
    return(list("chain_file"=mcmc_chain_file,"history_file"=infectionHistory_file,
                "step_scale"=steps))
}




