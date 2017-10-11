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
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")

    
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

    ## Initial indexing parameters
    no_recorded <- 1
    sampno <- 2
    par_i <- 1
    chain_index <- 1

    ## Table for storing infection histories
    historyTab <- matrix(NA, nrow=n_part*(save_block*25),ncol=n_strains+2)
