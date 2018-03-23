#' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. The algorithm then uses \code{\link{univ_proposal}} or \code{\link{mvr_proposal}} to explore parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename. This version of the algorithm is also designed to explore posterior densities for infection histories.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param mcmcPars named vector named vector with parameters for the MCMC procedure. Iterations (number of post adaptive iterations), popt (desired acceptance rate), popt_hist (desired acceptance rate for infection histories) opt_freq (after how many iterations do we adapt proposal), thin (save every n iterations), adaptive_period (number of iterations to adapt), save_block (number of post thinning iterations to save at a time), thin2 (infection history thinning), histSampleProb (proportion of inf histories to resample), switch_sample (resample inf histories every n iterations); burnin (number of iterations to run before any adapting), nInfs (number of infections to resample for each individual at each iteration), moveSizes (number of infection years/months to move when performing swap step), histProposal (which infection history proposal version to use, see \code{\link{describe_proposals}}, histOpt (if 1, performs adaptive infection history proposals. If 0, retains the starting infection history proposal parameters)
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param CREATE_POSTERIOR_FUNC pointer to posterior function used to calculate a likelihood
#' @param mvrPars leave NULL to use univariate proposals. Otherwise, a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @param PRIOR_FUNC user function of prior for model parameters. Should take parameter values only
#' @param version which version of the posterior function to use? See \code{\link{create_post_func}}
#' @param OPT_TUNING constant describing the amount of leeway when adapting the proposals steps to reach a desired acceptance rate (ie. does not change step size if within OPT_TUNING of the specified acceptance rate)
#' @param antigenicMap a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years 
#' @param ages data frame of ages and individual IDs for each participant, used to mask infection history proposals. Columns: age, individual. Can be left NULL
#' @param startInfHist infection history matrix to start MCMC at. Can be left NULL
#' @param ... other arguments to pass to CREATE_POSTERIOR_FUNC
#' @return a list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) a full file path at which the infection history chain is saved as a .csv file; 3) the last used covariance matrix; 4) the last used scale/step size (if multivariate proposals)
#' @export
run_MCMC_lambda<- function(parTab,data,mcmcPars,filename="test",OPT_TUNING=0.1,antigenicMap,ages,startInfHist,ver=1,mvrPars=NULL,block_weights=NULL,...){
    ## Extract MCMC parameters
    iterations <- mcmcPars["iterations"] # How many iterations to run after adaptive period
    popt <- mcmcPars["popt"] # Desired optimal acceptance rate
    popt_hist <- mcmcPars["popt_hist"]
    opt_freq<- mcmcPars["opt_freq"] # How often to adjust step size
    thin <- mcmcPars["thin"] # Save only every nth iterations for theta sampling
    adaptive_period<- mcmcPars["adaptive_period"] # How many iterations for adaptive period
    save_block <- mcmcPars["save_block"] # How many post-thinning iterations to store before saving to disk
    histTabThin <- mcmcPars["thin2"] # Save only every nth iterations for infection history sampling
    histSampleProb <- mcmcPars["histSampleProb"] # What proportion of infection histories to sample each step
    switch_sample <- mcmcPars["switch_sample"] # Resample infection histories every n iterations
    burnin <- mcmcPars["burnin"] # Run this many iterations before attempting adaptation. Idea is to reduce getting stuck in local maxima
    moveSize <- mcmcPars["moveSize"] # Number of infections to move/remove/add in each proposal step
    nInfs <- mcmcPars["nInfs"] # Number of infections to move/remove/add in each proposal step
    histOpt <- mcmcPars["histOpt"] # Should infection history proposal step be adaptive?
    n_years <- mcmcPars["nYears"]
    ## Extract parameter settings
    if(is.null(parTab$block)){
        unique_blocks <- c(0,1)
        blocks <- rep(1, nrow(parTab))
    } else {
        unique_blocks <- c(0,unique(parTab$block))
        if(!is.null(block_weights)) unique_blocks <- rep(unique_blocks, block_weights)
        print(unique_blocks)
        #unique_blocks <- c(0,2)
        blocks <- parTab$block
    }
    
    param_length <- nrow(parTab)
    unfixed_pars <- which(parTab$fixed == 0 & parTab$block == 1) # Indices of free parameters
    unfixed_par_length <- nrow(parTab[parTab$fixed== 0 & parTab$block == 1,]) # How many free parameters?
    current_pars <- parTab$values # Starting parameters
    par_names <- as.character(parTab$names) # Parameter names
    ## Parameter constraints
    lower_bounds <- parTab$lower_bound # Parameters cannot step below this
    upper_bounds <- parTab$upper_bound # Parameters cannot step above this
    steps <- parTab$steps # How far to step on unit scale to begin with?
    fixed <- parTab$fixed # Index which parameters are fixed

    lambda_indices <- which(parTab$names == "lambda")
    
    ## Pull out alpha and beta for beta binomial proposals
    if("alpha" %in% par_names & "beta" %in% par_names){
        alpha <- parTab[parTab$names == "alpha","values"]
        beta <- parTab[parTab$names == "beta","values"]
    }
    
    ## Arrays to store acceptance rates
    ## If univariate proposals, store vector of acceptances
    if(is.null(mvrPars)){
        tempaccepted <- tempiter <- integer(param_length)
        reset <- integer(param_length)
        reset[] <- 0
    } else {
        ## If multivariate proposals, aggregate to one acceptance rate.
        ## Also extract covariance matrix, scale of proposal steps, and how
        ## much weighting to give to previous covariance matrix upon adaptive update
        tempaccepted <- tempiter <- reset <- 0
        covMat <- mvrPars[[1]][unfixed_pars,unfixed_pars]
        covMat0 <- diag(unfixed_par_length)
        steps <- mvrPars[[2]]
        w <- mvrPars[[3]]
    }
    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste0(filename,"_chain.csv")
    infectionHistory_file <- paste0(filename,"_infectionHistories.csv")
    
###############
    ## Extract data parameters
##############
    strainIsolationTimes <- unique(antigenicMap$inf_years) # How many strains are we testing against and what time did they circulate
    samplingTimes <- unique(data$sample) # What are the range of sampling times?
    n_strain <- length(strainIsolationTimes) # How many strains could an individual see?
    n_indiv <- length(unique(data$individual)) # How many individuals in the data?
    individuals <- 1:n_indiv # Create vector of individuals

    ###################
    ## Housekeeping for infection history chain
    ###################
    histiter <- histaccepted <- histiter_add <- histaccepted_add <- histiter_move <- histaccepted_move <- histreset <- integer(n_indiv)
    nInfs_vec <- rep(nInfs, n_indiv) # How many infection history moves to make with each proposal
    moveSizes <- rep(moveSize, n_indiv) # How many years to move in smart proposal step
    
    ageMask <- create_age_mask(ages, strainIsolationTimes, n_indiv)
    n_alive <- sapply(strainIsolationTimes, function(x) nrow(ages[ages$DOB <= x,]) )
    ## Create posterior calculating function
    posterior_simp_full <- protect(create_post_func(parTab,data,antigenicMap,PRIOR_FUNC=NULL,4,ageMask,...))
    posterior_simp_data <- protect(create_post_func(parTab,data,antigenicMap,PRIOR_FUNC=NULL,1,ageMask,...))
    posterior_simp_lambda<- protect(create_post_func(parTab,data,antigenicMap,PRIOR_FUNC=NULL,3,ageMask,...))

######################
    ## Setup initial conditions
    infectionHistories = startInfHist
   
    ## Initial likelihood
    probabs_dat <- posterior_simp_data(current_pars,infectionHistories)
    probabs_lambda <- posterior_simp_lambda(current_pars,infectionHistories)
    #probabs_full <- posterior_simp_full(current_pars,infectionHistories)
    probabs_full <- probabs_dat + probabs_lambda
    probab_full <- sum(probabs_full)
    probab_dat <- sum(probabs_dat)
    probab_lambda <- sum(probabs_lambda)
    message(cat("Starting theta likelihood: ",probab_full,sep="\t"))
###############
    
####################
    ## PRE ALLOCATE MEMORY
####################
    ## Create empty chain to store every iteration for the adaptive period and burnin
    opt_chain <- matrix(nrow=burnin + adaptive_period,ncol=unfixed_par_length)

    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=param_length+2+2+1)

    ## Set up initial csv file
    chain_colnames <- c("sampno",par_names,"lnlike","prob_dat","prob_lambda","block")
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,current_pars,probab_full, probab_dat, probab_lambda,blocks[2])
    colnames(tmp_table) <- chain_colnames
    
    ## Write starting conditions to file
    data.table::fwrite(as.data.frame(tmp_table),file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)
    save_infHist_to_disk(infectionHistories, infectionHistory_file, 1, append=FALSE,colNames=TRUE)
    ## Initial indexing parameters
    no_recorded <- 1
    no_recorded_infHist <- 1
    sampno <- 2
    par_i <- 1
    chain_index <- 1
    last_index <- 1
    block_index <- 2
    block <- blocks[block_index]
#####################
    ## MCMC ALGORITHM
#####################
    
    for(i in 1:(iterations + adaptive_period + burnin)){
        block <- unique_blocks[block_index]
        if(i %% save_block == 0) message(cat("Current iteration: ", i, sep="\t"))
######################
        ## PROPOSALS
######################
        ## If updating theta
                                        #if(i %% switch_sample == 0){
        if(block == 1){
            ## If using univariate proposals
            ## For each parameter (Gibbs)
            ## If using univariate proposals
            if(is.null(mvrPars)) {
                j <- unfixed_pars[par_i]
                par_i <- par_i + 1
                if(par_i > unfixed_par_length) par_i <- 1
                proposal <- univ_proposal(current_pars, lower_bounds, upper_bounds, steps,j)
                tempiter[j] <- tempiter[j] + 1
            } else {
                ## NOTE
                ## MIGHT WANT TO USE ADAM'S PROPOSAL FUNCTION
                tmp <- mvr_proposal(current_pars, unfixed_pars, steps*covMat,
                                    covMat0, TRUE, beta=0.05,lower=lower_bounds, upper=upper_bounds)
                proposal <- tmp
                tempiter <- tempiter + 1
            }
            ## Calculate new likelihood for these parameters
            new_probabs_dat <- posterior_simp_data(proposal,infectionHistories)
            new_probabs_lambda<- posterior_simp_lambda(proposal,infectionHistories)
            new_probab_dat <- sum(new_probabs_dat)
            new_probab_lambda <- sum(new_probabs_lambda)
            new_probabs_full <- new_probabs_dat + new_probabs_lambda
            new_probab_full <- new_probab_dat + new_probab_lambda
            
            ## Otherwise, resample infection history
        } else if(block == 2){
            years <- 1:ncol(infectionHistories)
            years <- sample(years, n_years)
            js <- lambda_indices[years]
            tmp <- lambda_proposal(current_pars, infectionHistories, years, js, alpha, beta, n_alive)
            proposal <- tmp[[1]]
            ratio <- tmp[[2]]
            ## Calculate new likelihood for these parameters
            new_probabs_lambda <- posterior_simp_lambda(proposal,infectionHistories)            
            new_probab_lambda <- sum(new_probabs_lambda)
        } else {
            ## Choose a random subset of individuals to update
            indivSubSample <- sample(1:n_indiv, ceiling(histSampleProb*n_indiv))
            randNs <- runif(length(indivSubSample))
            ## Which infection history proposal to use?
            
            newInfectionHistories <- infection_history_symmetric(infectionHistories, indivSubSample, ageMask, moveSizes, nInfs_vec, randNs)
            
            ## The proposals are either a swap step or an add/remove step. Need to track which type was used for which individual,
            ## as we adapt the `step size` for these two update steps independently
            move <- which(randNs > 1/2)
            add <- which(randNs < 1/2)
            
            histiter_add[indivSubSample[add]]<- histiter_add[indivSubSample[add]] + 1
            histiter_move[indivSubSample[move]]<- histiter_move[indivSubSample[move]] + 1
            
            ## Calculate new likelihood with these infection histories
            new_probabs_dat<- posterior_simp_data(current_pars, newInfectionHistories)
            new_probab_dat<- sum(new_probabs_dat)

            new_probabs_lambda<- posterior_simp_lambda(current_pars, newInfectionHistories)
            new_probab_lambda<- sum(new_probabs_lambda)
            
            new_probabs_full <- new_probabs_dat + new_probabs_lambda
            new_probab_full <- new_probab_dat + new_probab_lambda
            
            histiter[indivSubSample]<- histiter[indivSubSample] + 1
        }

#############################
        ## METROPOLIS HASTINGS STEP
#############################
        ## Check that all proposed parameters are in allowable range
        ## Skip if any parameters are outside of the allowable range
        if(block == 1){
            log_prob <- new_probab_dat - probab_dat
            log_prob <- min(log_prob, 0)
            if((is.finite(log_prob) && log(runif(1)) < log_prob)){
                if(!any(proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                        proposal[unfixed_pars] > upper_bounds[unfixed_pars])){
                    ## Accept with probability 1 if better, or proportional to
                    ## difference if not
                    current_pars <- proposal
                    ## Store acceptances
                    tempaccepted[j] <- tempaccepted[j] + 1
                    probabs_dat <- new_probabs_dat
                    probab_dat <- new_probab_dat
                    probabs_lambda <- new_probabs_lambda
                    probab_lambda <- new_probab_lambda
                    probabs_full <- new_probabs_dat + new_probabs_lambda
                    probab_full <- new_probab_dat + new_probab_lambda
                }
            }
        } else if(block == 2){
            #log_prob <- new_probab_lambda - probab_lambda + ratio
                                        #message(cat("ratio: ", ratio,sep="\t"))
            #log_prob <- min(log_prob, 0)
#            print(log_prob)
            #if((is.finite(log_prob) && log(runif(1)) < log_prob)){
                current_pars <- proposal
                ## Store acceptances
                                        #probabs_dat <- probabs_dat
                                        #probab_dat <- probab_dat
                probabs_lambda <- new_probabs_lambda
                probab_lambda <- new_probab_lambda
                probabs_full <- probabs_dat + new_probabs_lambda
                probab_full <- probab_dat + new_probab_lambda
            #}
        } else {
            
            ## MH step for each individual
            log_probs <- (new_probabs_full[indivSubSample] - probabs_full[indivSubSample])#  + log(acceptance[indivSubSample])
            log_probs[log_probs > 0] <- 0
            x <- which(log(runif(length(indivSubSample))) < log_probs)
            changeI <- indivSubSample[x]
            infectionHistories[changeI,] <- newInfectionHistories[changeI,]
            probabs_full[changeI] <- new_probabs_full[changeI]
            probabs_dat[changeI] <- new_probabs_dat[changeI]
            probabs_lambda[changeI] <- new_probabs_lambda[changeI]
            probab_full<- sum(probabs_full)
            probab_dat <- sum(probabs_dat)
            probab_lambda <- sum(probabs_lambda)
            
            
            ## Record acceptances for each add or move step
            add <- intersect(add, changeI)
            move <- intersect(move, changeI)
            histaccepted_add[indivSubSample[add]] <- histaccepted_add[indivSubSample[add]] + 1
            histaccepted_move[indivSubSample[move]] <- histaccepted_move[indivSubSample[move]] + 1
            histaccepted[changeI] <- histaccepted[changeI] + 1
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
            save_chain[no_recorded,2:(ncol(save_chain)-1 - 3)] <- current_pars
            save_chain[no_recorded,(ncol(save_chain)-3):ncol(save_chain)] <- c(probab_full, probab_dat, probab_lambda,block)
            no_recorded <- no_recorded + 1
        }

        ## Save infection histories
        if(i %% histTabThin == 0){
           save_infHist_to_disk(infectionHistories, infectionHistory_file, sampno)
        }
        
##############################
        ## ADAPTIVE PERIOD
##############################
        ## If within adaptive period, need to do some adapting!
        if(i > (adaptive_period + burnin) & i %% opt_freq == 0){
            pcur <- tempaccepted/tempiter ## get current acceptance rate
            message(cat("Pcur: ", pcur,sep="\t"))
            message(cat("Step sizes: ", steps,sep="\t"))
            tempaccepted <- tempiter <- reset

            ## Have a look at the acceptance rates for infection histories
            pcurHist <- histaccepted/histiter ## Overall
            pcurHist_add <- histaccepted_add/histiter_add ## For adding
            pcurHist_move <- histaccepted_move/histiter_move ## For moving
            histiter <- histaccepted <- histaccepted_add <- histaccepted_move <- histiter_add <- histiter_move <- histreset
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
                message(cat("Hist acceptance add: ", pcurHist_add, cat="\t"))
                message(cat("Hist acceptance move: ", pcurHist_move, cat="\t"))
                message(cat("Mean hist acceptance: ", mean(pcurHist),cat="\t"))
                histiter <- histaccepted <- histaccepted_add <- histaccepted_move <- histiter_add <- histiter_move <- histreset

                ## If adaptive infection history proposal
                if(histOpt == 1){
                    nInfs_vec[which(pcurHist_add < popt_hist*(1-OPT_TUNING))] <- nInfs_vec[which(pcurHist_add< popt_hist*(1-OPT_TUNING))] - 1
                    nInfs_vec[which(pcurHist >= popt_hist*(1+OPT_TUNING))] <- nInfs_vec[which(pcurHist >= popt_hist*(1+OPT_TUNING))] +1
                    nInfs_vec[nInfs_vec < 1] <- 1

                    moveSizes[which(pcurHist < popt_hist*(1-OPT_TUNING))] <- moveSizes[which(pcurHist < popt_hist*(1-OPT_TUNING))] - 1
                    moveSizes[which(pcurHist >= popt_hist*(1+OPT_TUNING))] <- moveSizes[which(pcurHist >= popt_hist*(1+OPT_TUNING))] +1
                    moveSizes[moveSizes < 1] <- 1

                    for(ii in seq_along(nInfs_vec)){
                        moveSizes[ii] <- min(moveSizes[ii], n_strain - ageMask[ii])
                        nInfs_vec[ii] <- min(nInfs_vec[ii],n_strain - ageMask[ii])
                    }
                }
                message(cat("nInfs: ", nInfs_vec, sep="\t"))
                message(cat("Move sizes: ", moveSizes, sep="\t"))
                message(cat("Pcur: ", pcur,sep="\t"))
                message(cat("Step sizes: ", steps,sep="\t"))
                tempaccepted <- tempiter <- reset
                message(cat("Hist sample prob: ", histSampleProb))
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
        block_index <- block_index + 1
        if(block_index > length(unique_blocks)) block_index <- 1
        sampno <- sampno + 1
    }

    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        data.table::fwrite(as.data.frame(save_chain[1:(no_recorded-1),]),file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }
   
    return(list("chain_file"=mcmc_chain_file,"history_file"=infectionHistory_file,"step_scale"=steps))
}




