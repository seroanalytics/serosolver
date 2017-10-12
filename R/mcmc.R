#' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{univ_proposal}} or \code{\link{mvr_proposal}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param data the data frame of data to be fitted
#' @param mcmcPars named vector named vector with parameters for the MCMC procedure. Iterations, popt, opt_freq, thin, burnin, adaptive_period and save_block.
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param CREATE_POSTERIOR_FUNC pointer to posterior function used to calculate a likelihood
#' @param mvrPars a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @param PRIOR_FUNC user function of prior for model parameters. Should take values, names and local from param_table
#' @param OPT_TUNING constant used to indicate what proportion of the adaptive period should be used to build the covariance matrix, if needed
#' @param antigenicMap xxxx
#' @param samplingInformation xxxx
#' @param ages xxxx
#' @return a list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) a full file path at which the infection history chain is saved as a .csv file; 3) the last used covarianec matrix; 4) the last used scale size; 5) the last used step size (if univariate proposals)
#' @export
run_MCMC <- function(parTab,
                     data,
                     mcmcPars,
                     filename,
                     CREATE_POSTERIOR_FUNC,
                     mvrPars,
                     PRIOR_FUNC,
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

    ## Extract parameter settings
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

    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste0(filename,"_chain.csv")
    infectionHistory_file <- paste0(filename,"_infectionHistories.csv")
    
###############
    ## Extract data parameters
##############
    strainIsolationTimes <- unique(data$strain)
    samplingTimes <- unique(data$sample)
    n_strains <- length(strainIsolationTimes)
    n_indiv <- length(unique(data$individual))
    n_groups <- length(unique(data$group))
    individuals <- 1:n_indiv
    groups <- 1:n_groups
    
    ## Create posterior calculating function
    posterior_simp <- protect(CREATE_POSTERIOR_FUNC(parTab,data,antigenicMap,
                                                    sampingInformation,
                                                    PRIOR_FUNC,...))
   
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

######################
    ## Setup initial conditions
    infectionHistories <- setup_infection_histories(data, strainIsolationTimes, ageMask)
    
    ## Initial likelihood
    #### NOTE
    #### IT MIGHT BE A BIT TOO SLOW PASSING INFECTION HISTORIES AS A MATRIX EACH ITERATION
    probab <- posterior_simp(current_pars,infectionHistories)
###############
    
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
    historyTab <- emptyHistTab <- matrix(NA, nrow=n_indiv*(save_block*n_indiv),ncol=n_strains+2)
    tmp_table <- matrix(NA, nrow=n_indiv,ncol=n_strains+2)
    tmp_table[1:nindiv,1:n_strain] <- infectionHistories
    tmp_table[1:nindiv,n_strain+1] <- individuals
    tmp_table[1:nindiv,n_strain+2] <- 1
    colnames(tmp_table) <- c(as.character(strainIsolationTimes),"individual","sampno")
    ## Write starting infectoin histories
    write.table(tmp_table, infectionHistory_file, row.names=FALSE, col.names=FALSE, sep=",",append=FALSE)

    ## Initial indexing parameters
    no_recorded <- 1
    no_recorded_infhist <- 1
    sampno <- 2
    par_i <- 1
    chain_index <- 1

#####################
    ## MCMC ALGORITHM
#####################
    for(i in 1:(iterations + adaptive_period)){
        if(i %% save_block == 0) message(cat("Current iteration: ", i, sep="\t"))
######################
        ## PROPOSALS
######################
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
                ## NOTE
                ## MIGHT WANT TO USE ADAM'S PROPOSAL FUNCTION
                proposal <- mvr_proposal(current_pars, unfixed_pars, scale*covMat)
                tempiter <- tempiter + 1
            }
            ## Calculate new likelihood for these parameters
            new_probab <- posterior_simp(proposal, infectionHistories)
            ## Otherwise, resample infection history
        } else {
            indivSubSample <- sample(1:n_indiv, ceiling(histSampleProb*n_indiv))
            newInfectionHistory <- infection_history_proposal(infectionHistories, indivSubSample,
                                                              strainIsolationTimes, ageMask)
            ## Calculate new likelihood with these infection histories
            new_probab <- posterior_simp(current_pars, newInfectionHistories)
        }
#########################
        ## We could add a function pointer that does the same job as "pmask" - this way it is a bit more
        ## generalised. Shouldn't need to though, thanks to the function pointer idea
#########################

        #############################
        ## METROPOLIS HASTINGS STEP
        #############################
        ## Check that all proposed parameters are in allowable range
        ## Skip if any parameters are outside of the allowable range
        if(!any(proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                proposal[unfixed_pars] > upper_bounds[unfixed_pars])){
            
            log_prob <- min(new_probab-probab,0)
            ## Accept with probability 1 if better, or proportional to
            ## difference if not
            if(is.finite(log_prob) && log(runif(1)) < log_prob){
                if(i %% switch_sample == 0){
                    current_pars <- proposal
                    ## Store acceptances
                    if(is.null(mvrPars)) tempaccepted[j] <- tempaccepted[j] + 1
                    else tempaccepted <- tempaccepted + 1
                } else {
                    infectionHistories <- newInfectionHistories
                }
                probab <- new_probab
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
            save_chain[no_recorded,2:(ncol(save_chain)-1)] <- current_pars
            save_chain[no_recorded,ncol(save_chain)] <- probab
            no_recorded <- no_recorded + 1
        }

        ## Save infection histories
        if(i %% histTabThin == 0){
            historyTab[no_recorded_infHist:(no_recorded_infHist + n_indiv-1),1:n_strain] <- infectionHistories
            historyTab[no_recorded_infHist:(no_recorded_infHist + n_indiv-1),n_strain+1] <- individuals
            historyTab[no_recorded_infHist:(no_recorded_infHist + n_indiv-1),n_strain+2] <- i
            no_recorded_infHist <- no_recorded_infHist + n_indiv
        }
        
##############################
        ## ADAPTIVE PERIOD
##############################
        ## If within adaptive period, need to do some adapting!
        if(i <= adaptive_period){
            ## Current acceptance rate
            pcur <- tempaccepted/tempiter
            ## Save each step
            opt_chain[chain_index,] <- current_pars[unfixed_pars]
            ## If in an adaptive step
            if(chain_index %% opt_freq == 0){
                ## If using univariate proposals
                if(is.null(mvrPars)){
                    ## For each non fixed parameter, scale the step size
                    for(x in unfixed_pars) steps[x] <- scaletuning(steps[x],popt,pcur[x])
                    message(cat("Pcur: ", pcur[unfixed_pars],sep="\t"))
                    message(cat("Step sizes: ", steps[unfixed_pars],sep="\t"))
                    tempaccepted <- tempiter <- reset
                } else {       ## If using multivariate proposals
                    if(chain_index > OPT_TUNING*adaptive_period & chain_index < (0.8*adaptive_period)){
                        oldCovMat <- covMat
                        ## Creates a new covariance matrix, but weights it with the old one
                        covMat <- cov(opt_chain[1:chain_index,])
                        covMat <- w*covMat + (1-w)*oldCovMat
                    }
                    ## Scale tuning for last 20% of the adpative period
                    if(chain_index > (0.8)*adaptive_period){
                        scale <- scaletuning(scale, popt,pcur)
                    }
                    tempiter <- tempaccepted <- 0
                    message(cat("Pcur: ", pcur,sep="\t"))
                    message(cat("Scale: ", scale,sep="\t"))
                }
            }
            chain_index <- chain_index + 1
        }
        #######################
        ## HOUSEKEEPING
#######################
        if(no_recorded == save_block){
            write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,
                        col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            save_chain <- empty_save_chain
            no_recorded <- 1
        }
        if(no_recorded_infHist/n_indiv == save_block){
            write.table(historyTab[1:(no_recorded_infHist-n_indiv),], file=infectionHistory_file,
                        col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            historyTab <- emptyHistoryTab
            no_recorded_infHist <- 1
        }
        sampno <- sampno + 1
    }

    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }
    if(no_recorded_infHist > 2){
        write.table(historyTab, file=infectionHistory_file,
                    col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
        historyTab <- emptyHistoryTab
        no_recorded_infHist <- 1
    }

    if(is.null(mvrPars)){
        covMat <- NULL
        scale <- NULL
    } else {
        steps <- NULL
    }
    return(list("chain_file"=mcmc_chain_file,"history_file"=infectionHistory_file,
                "covMat"=covMat,"scale"=scale, "steps"=steps))
}




