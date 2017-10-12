run_MCMC <- function(parTab,
                     data,
                     mcmcPars,
                     filename,
                     CREATE_POSTERIOR_FUNC,
                     mvrPars,
                     PRIOR_FUNC
                     OPT_TUNING,
                     antigenicMap,
                     samplingInformation,
                     ages=NULL,
                     ...){
    ## Extract MCMC parameters
    iterations <- mcmcPars["iterations"]
    popt <- mcmcPars["popt"]
    opt_freq<- mcmcPars["opt_freq"]
    thin <- mcmcPars["thin"]
    adaptive_period<- mcmcPars["adaptive_period"]
    save_block <- mcmcPars["save_block"]
    histTabThin <- mcmcPars["thin2"]
    histSampleProb <- mcmcPars["histSampleProb"]
    switch_sample <- mcmcPars["switch_sample"]
    
    param_length <- nrow(parTab)
    unfixed_pars <- which(parTab$fixed == 0)
    unfixed_par_length <- nrow(parTab[parTab$fixed== 0,])
    current_pars <- parTab$values
    par_names <- as.character(parTab$names)

    ## Parameter constraints
    lower_bounds <- parTab$lower_bound
    upper_bounds <- parTab$upper_bound
    steps <- parTab$steps
    fixed <- parTab$fixed
    
    ## Arrays to store acceptance rates
    ## If univariate proposals
    if(is.null(mvrPars)){
        tempaccepted <- tempiter <- integer(param_length)
        reset <- integer(param_length)
        reset[] <- 0
    } else { # If multivariate proposals
        tempaccepted <- tempiter <- 0
        covMat <- mvrPars[[1]][unfixed_pars,unfixed_pars]
        scale <- mvrPars[[2]]
        w <- mvrPars[[3]]
    }

    ## Create posterior calculating function
    posterior_simp <- protect(CREATE_POSTERIOR_FUNC(parTab,data, 
                                                    PRIOR_FUNC,...))
    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste0(filename,"_chain.csv")
    infectionHistory_file <- paste0(filename,"_infectionHistories.csv")
    
###############
    ## Read in antigenic map and extract data parameters
##############
    antigenicMapMelted <- c(outputdmatrix.fromcoord(antigenicMap))
    strainIsolationTimes <- unique(dat$strain)
    samplingTimes <- unique(dat$sample)
    n_strains <- length(strainIsolationTimes)
    
###############
    ## Create age mask
    ## Note that ages for all groups must be from same reference point
###############
    if(!is.null(ages)){
        tref <- max(max(strainIsolationTimes),max(samplingTimes))
        ages$DOB <- tref - ages$age
        ageMask <- sapply(ages$DOB, function(x){
            if(is.na(x)){
                1
            } else {
                which(as.numeric(x <= strainIsolationTimes) > 0)[1]
            }
        })
    } else {
        ageMask <- rep(1, n_indiv)
    }

####################
    ## PRE ALLOCATE MEMORY
####################
    ## Create empty chain to store every iteration for the adaptive period
    opt_chain <- matrix(nrow=adaptive_period,ncol=unfixed_par_length)

    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=param_length+2+misc_length)

    ## Set up initial csv file
    chain_colnames <- c("sampno",par_names,misc_colnames,"lnlike")
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,current_pars,misc,probab)
    colnames(tmp_table) <- chain_colnames
    ## Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    ## Table for storing infection histories
    historyTab <- matrix(NA, nrow=n_part*(save_block*25),ncol=n_strains+2)
    historyTab[1:25,1:n_strain] <- infectionHistories
    historyTab[,n_strain+1] <- 1:n_indiv
    historyTab[,n_strain+2] <- 1
    write.table(historyTab, infectionHistory_file, row.names=FALSE, col.names=FALSE, sep=",",append=FALSE)


    ## Initial indexing parameters
    no_recorded <- 1
    sampno <- 2
    par_i <- 1
    chain_index <- 1

    #####################
    ## MCMC ALGORITHM
#####################
    
    for(i in 1:(iterations + adaptive_period)){
        ## If updating theta
        if(i %% switch_sample == 0){
            ## If using univariate proposals
            if(is.null(mvrPars)) {
                ## For each parameter (Gibbs)
                j <- unfixed_pars[par_i]
                par_i <- par_i + 1
                if(par_i > unfixed_par_length) par_i <- 1
                proposal <- univ_proposal(current_pars, lower_bounds, upper_bounds, steps,j)
                tempiter[j] <- tempiter[j] + 1
                ## If using multivariate proposals
            } else {
                proposal <- mvr_proposal(current_pars, unfixed_pars, scale*covMat)
                tempiter <- tempiter + 1
            }
            
            ## Otherwise, resample infection history
        } else {
            indivSubSample <- sample(1:n_indiv, ceiling(histSampleProb*n_indiv))
            newInfectionHistory <- infection_history_proposal(infectionHistory, indivSubSample,
                                                              strainIsolationTimes, ageMask)
        }
        antigenicMapLong <- 1 - current_pars["sigma"]*antigenicMapMelted
        antigenicMapLong[antigenicMapLong < 0] <- 0
        antigenicMapShort <- 1 - current_pars["sigma2"]*antigenicMapMelted
        antigenicMapShort[antigenicMapShort < 0] <- 0
        

        ## Propose new parameters and calculate posterior
        ## Check that all proposed parameters are in allowable range
        ## Skip if any parameters are outside of the allowable range
        if(!any(
                proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                proposal[unfixed_pars] > upper_bounds[unfixed_pars]
            )
           ){
            ## Calculate new likelihood and find difference to old likelihood
            new_probab <- posterior_simp(proposal)
            log_prob <- min(new_probab-probab,0)
            
            ## Accept with probability 1 if better, or proportional to
            ## difference if not
            if(is.finite(log_prob) && log(runif(1)) < log_prob){
                current_pars <- proposal
                probab <- new_probab
                
                ## Store acceptances
                if(is.null(mvrPars)){
                    tempaccepted[j] <- tempaccepted[j] + 1
                } else {
                    tempaccepted <- tempaccepted + 1
                }
            }
        }









        
    }
    
    
    
