run_MCMC2 <- function(parTab,data,mcmcPars,filename="test",CREATE_POSTERIOR_FUNC,PRIOR_FUNC=NULL,version=1,
                      OPT_TUNING=0.1,antigenicMap,ages=NULL,startInfHist=NULL,...){
    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste0(filename,"_chain.csv")
    infectionHistory_file <- paste0(filename,"_infectionHistories.csv")
    
################################
    ## Extract MCMC parameters
################################
    iterations <- mcmcPars["iterations"] # How many iterations to run after adaptive period
    opt_freq<- mcmcPars["opt_freq"] # How often to adjust step size
    thin <- mcmcPars["thin"] # Save only every nth iterations for theta sampling
    adaptive_period<- mcmcPars["adaptive_period"] # How many iterations for adaptive period
    save_block <- mcmcPars["save_block"] # How many post-thinning iterations to store before saving to disk
    histTabThin <- mcmcPars["thin2"] # Save only every nth iterations for infection history sampling
    burnin <- mcmcPars["burnin"] # Run this many iterations before attempting adaptation. Idea is to reduce getting stuck in local maxima
    histProposal <- mcmcPars["histProposal"] # Which infection history proposal version?
    
################################
    ## Extract parameter settings
################################
    current_pars <- parTab$values # Starting parameters
    param_length <- nrow(parTab)
    unfixed_pars <- which(parTab$fixed == 0) # Indices of free parameters
    unfixed_par_length <- nrow(parTab[parTab$fixed== 0,]) # How many free parameters?
    par_names <- as.character(parTab$names) # Parameter names
    blocks <- parTab$block # Updating parameters in blocks
    unique_blocks <- unique(blocks) # Need list of blocks to cycle through
    inf_block <- max(unique_blocks) + 1 # Add an extra block for updating infection histories
    ## If we want to weight sampling of blocks, just add more of that block to "unique_blocks"
    unique_blocks <- c(unique_blocks, inf_block)
    n_blocks <- length(unique_blocks)
    
###############
    ## Extract data parameters
##############
    strainIsolationTimes <- unique(antigenicMap$inf_years) # How many strains are we testing against and what time did they circulate
    samplingTimes <- unique(data$sample) # What are the range of sampling times?
    n_strain <- length(strainIsolationTimes) # How many strains could an individual see?
    n_indiv <- length(unique(data$individual)) # How many individuals in the data?
    individuals <- 1:n_indiv # Create vector of individuals

###############
    ## Create age mask
    ## -----------------------
    ## Note that ages for all groups must be from same reference point
    ## -----------------------
###############
    if(!is.null(ages)){
        ageMask <- create_age_mask(ages, strainIsolationTimes, n_indiv)
    } else {
        ageMask <- rep(1, n_indiv)
    }

    ## Setup initial conditions
    infectionHistories = startInfHist
    if(is.null(startInfHist)){
        infectionHistories <- setup_infection_histories_new(data, ages, strainIsolationTimes,
                                                            space=5,titre_cutoff=3)
    }
     ## Create posterior calculating function
    posterior_simp <- protect(CREATE_POSTERIOR_FUNC(parTab,data,antigenicMap,PRIOR_FUNC,version, ageMask,...))
    probabs <- posterior_simp(current_pars,infectionHistories)
    probab <- sum(probabs)
    message(cat("Starting theta likelihood: ",probab,sep="\t"))

    ## Set up unique proposal bundle, proposal function and proposal steps for each
    ## parameter block
    proposal_bundle <- vector(mode="list",n_blocks)
    proposal_funcs <- vector(mode="list",n_blocks)
    acceptance_funcs <- vector(mode="list",n_blocks)
    for(block in unique_blocks){
        if(block != inf_block){
            proposal_ver <- unique(parTab[parTab$block==block,"proposal"])
            posterior_ver <- unique(parTab[parTab$block==block,"proposal"])
            proposal_bundle[[block]] <- create_proposal_bundle(parTab, proposal_ver, posterior_ver)
            proposal_funcs[[block]] <- create_proposal_funcs(parTab, proposal_ver, posterior_ver)
            acceptance_funcs[[block]] <- create_acceptance_funcs(parTab, proposal_ver, posterior_ver)
        } else {
            proposal_bundle[[block]] <- create_proposal_bundle_infHist(parTab, infectionHistories, mcmcPars)
            proposal_funcs[[block]] <- create_proposal_func_infHist(parTab, infectionHistories, mcmcPars)
            acceptance_funcs[[block]] <- create_acceptance_func_infHist(parTab, infectionHistories)                
        }
    }
    
####################
    ## PRE ALLOCATE MEMORY
####################
    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=param_length+2)

    ## Set up initial csv file
    chain_colnames <- c("sampno",par_names,"lnlike")
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,current_pars,probab)
    colnames(tmp_table) <- chain_colnames
    
    ## Write starting conditions to file
    data.table::fwrite(as.data.frame(tmp_table),file=mcmc_chain_file,row.names=FALSE,
                       col.names=TRUE,sep=",",append=FALSE)
    save_infHist_to_disk(infectionHistories, infectionHistory_file, 1, append=FALSE,colNames=TRUE)
    
    ## Initial indexing parameters
    no_recorded <- 1
    sampno <- 2

    ## First block
    block_index <- 1
    cur_block <- unique_blocks[block_index]

#####################
    ## MCMC ALGORITHM
#####################
    for(i in 1:(iterations + adaptive_period + burnin)){
        if(i %% save_block == 0) message(cat("Current iteration: ", i, sep="\t"))
######################
        ## PROPOSALS AND METROPOLIS HASTINGS STEP
######################
        if(cur_block != inf_block){
            ## Need to update number of iterations
            proposal_bundle[[cur_block]] <- proposal_funcs[[cur_block]](current_pars, infectionHistories, proposal_bundle[[cur_block]])
            proposal <- acceptance_funcs[[cur_block]](proposal_bundle, infectionHistories, probab, probabs)
            current_pars <- proposal[["pars"]]
            probab <- proposal[["new_probab"]]
            probabs <- proposal[["new_probabs"]]
        } else {
            proposal_bundle[[cur_block]] <- proposal_funcs[[cur_block]](current_pars, infectionHistories, proposal_bundle[[cur_block]])
            proposal <- acceptance_funcs[[cur_block]](proposal_bundle, infectionHistories, probab, probabs)
            infectionHistories <- proposal[["newInfectionHistories"]]
            probab <- proposal[["new_probab"]]
            probabs <- proposal[["new_probabs"]]
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
            save_chain[no_recorded,2:(ncol(save_chain)-1)] <- current_pars
            save_chain[no_recorded,ncol(save_chain)] <- probab
            no_recorded <- no_recorded + 1
        }
        ## Save infection histories
        if(i %% histTabThin == 0){
            save_infHist_to_disk(infectionHistories, infectionHistory_file, sampno)
        }
        
##############################
        ## ADAPTIVE PERIOD
##############################
        ## If after adaptive period, print status of each proposal bundle
        if(i > (adaptive_period + burnin) & i %% opt_freq == 0) monitor_status(proposal_bundle)
        if(i > burnin & i <= (adaptive_period + burnin) & i %% opt_freq == 0) proposal_bundle <- adaptive_update(proposal_bundle)
        
#######################
        ## HOUSEKEEPING
#######################
        if(no_recorded == save_block){
            data.table::fwrite(as.data.frame(save_chain[1:(no_recorded-1),]),file=mcmc_chain_file,
                               col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            save_chain <- empty_save_chain
            no_recorded <- 1
        }
        sampno <- sampno + 1
        block_index <- block_index + 1
        if(block_index > n_blocks) block_index <- 1
        cur_block <- unique_blocks[block_index]
    }
    
    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        data.table::fwrite(as.data.frame(save_chain[1:(no_recorded-1),]),file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }

    return(list("chain_file"=mcmc_chain_file,"history_file"=infectionHistory_file,"step_scale"=steps))
}
