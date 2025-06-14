#' Check infection history matrix
#'
#' Checks that the infection history matrix is allowable given the birth dates and sampling times of the data
#' @param antibody_data the data frame of titre data
#' @param possible_exposure_times vector of times at which individuals could be exposed e.g., `seq(1968,2015,by=1)`
#' @param verbose if TRUE, prints warning messages
#' @param inf_hist the infection history matrix, with nrows = number indivs and ncols = length(possible_exposure_times)
#' @return a single boolean value, FALSE if the check passes, otherwise throws an error
#' @family check_inputs
#' data(example_antibody_data)
#' data(example_inf_hist)
#' times <- 1968:2015
#' check_inf_hist(example_antibody_data, times, example_inf_hist)
#' @export
check_inf_hist <- function(antibody_data, possible_exposure_times, inf_hist,verbose=FALSE){
    DOBs <- get_DOBs(antibody_data)
    age_mask <- create_age_mask(DOBs[,2],possible_exposure_times)
    sample_mask <- create_sample_mask(antibody_data, possible_exposure_times)
    before_born <- logical(length(age_mask))
    after_sample <- logical(length(sample_mask))
    res <- logical(length(age_mask))
    for(i in seq_along(age_mask)){
        if(sum(inf_hist[i,] != 0)){
            first_inf <- min(which(inf_hist[i,] == 1))
            last_inf <- max(which(inf_hist[i,] == 1))
            
            before_born[i] <- first_inf < age_mask[i]
            after_sample[i] <- last_inf > sample_mask[i]
            res[i] <- before_born[i] | after_sample[i]
        } else {
            res[i] <- FALSE
        }
    }
    if (any(res)) {
      if(verbose){
        message(cat("Which infections before birth: ", which(before_born), "\n"))
        message(cat("Which infections after last sample: ", which(after_sample), "\n"))
      }
        stop("Error in inf hist - infections occuring before individuals are born of after their latest sample\n")
    }
    return(any(res))
}


#' Check par_tab for simulate_data
#'
#' Checks the entries of par_tab used in simulate_data
#' @param par_tab the parameter table controlling information such as bounds, initial values etc
#' @param mcmc logical, if TRUE then checks are performed for the MCMC algorithm. Use FALSE when simulating data
#' @param version which version of the posterior function is being used? See \code{\link{create_posterior_func}}
#' @param possible_exposure_times optional vector of possible exposure times
#' @param verbose if TRUE, prints warning messages
#' @return the same par_tab object with corrections if needed
#' @family check_inputs
#' @examples
#' data(example_par_tab)
#' check_par_tab(example_par_tab, FALSE, version=1)
#' @export
check_par_tab <- function(par_tab, mcmc = FALSE, version = NULL,possible_exposure_times=NULL, verbose=FALSE) {
    ## Checks that should happen in simulate_data and serosolver
    essential_names <- c("names","values","fixed","lower_bound","upper_bound","lower_start","upper_start","par_type")
    if (!all(essential_names %in% colnames(par_tab))) {
        message(paste(c("Some column names missing from par_tab: ", setdiff(essential_names,colnames(par_tab)),"\n"),collapse=" "))
     
        if(!("type" %in% colnames(par_tab))){
          if(verbose) message(cat("Adding \"par_type\" to par_tab variables.\n"))
          par_tab$par_type <- 1
        }
    }
    pars <- par_tab$values
    names(pars) <- par_tab$names

    explicit_phi <- "phi" %in% names(pars)

    ## Additional checks that should happen in serosolver
    if (mcmc == TRUE) {
      
      if(!("steps" %in% colnames(par_tab))){
        if(verbose) message(cat("Adding \"steps\" to par_tab variables.\n"))
        par_tab$steps <- 0.1
      }
      
        phi_indices <- which(par_tab$names == "phi")
        no_phi <- length(phi_indices)
        explicit_phi <- (no_phi > 0)

        ## Check that all optional parameters are fixed, if not, fix them
        op_pars <- par_tab[which(par_tab$par_type == 0), ]
        if (all(op_pars$fixed == 1) == FALSE) stop("All optional parameters must be fixed")

        if (version == 1) {
            ## Check that the correct number of phis are present
            if(!explicit_phi){
                if(verbose) message(cat("Prior version 1 is specified, but there are no entries for phi in par_tab. Note that there must be one phi entry per entry in possible_exposure_times.\n"))
            }
          
            if(!is.null(possible_exposure_times) & no_phi < length(possible_exposure_times)){
              if(verbose) message(cat("Padding par_tab with phis, as number of phi parameters does not match possible_exposure_times.\n")) ## Should we remove them?
              par_tab_phi <- data.frame(names="phi",values=0.1,fixed=0,lower_bound=0,upper_bound=1,lower_start=0,upper_start=1,par_type=2,steps=0.1,stratification=NA,biomarker_group=NA)
              n_missing <- length(possible_exposure_times) - no_phi
              for(i in 1:n_missing){
                par_tab <- rbind(par_tab, par_tab_phi)
              }
            }
        }
        
        if (all(par_tab$fixed == 1)) {
            stop("All parameters are set to be fixed. The MCMC procedure will not work. Please set at least the entry for error to be fixed <- 0.\n")
        }

        
        if (version %in% c(2, 3, 4)) {
            if (explicit_phi){
              if(verbose) message(cat("Phis are not required for versions 2, 3 or 4 but par_tab contains phi(s). Removing all entries named phi.\n")) ## Should we remove them?
              par_tab <- par_tab[par_tab$names != "phi",]
            }
        }
        ## Check bounds are equal to starting bounds
        if (any(par_tab$upper_start > par_tab$upper_bound) | any(par_tab$lower_start < par_tab$lower_bound)) {
            warning("lower_start and upper_start are not equal to the starting lower_bound and upper_bound. If par_tab was used to create starting values, starting values may be out of bounds.\n ")
        }
    }
    ## Check that alpha and beta there for beta distribution
    ## If there, Pull out alpha and beta for beta binomial proposals
    if (!("infection_model_prior_shape1" %in% par_tab$names) | !("infection_model_prior_shape2" %in% par_tab$names)) {
        stop("par_tab must have entries for `infection_model_prior_shape1` and `infection_model_prior_shape2` for infection history prior.\n")
    }
    par_tab
}

#' Checks the entries of data used in serosolver
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param verbose if TRUE, prints warning messages
#' @return the same data object with corrections if needed
#' @family check_inputs
#' @examples
#' data(example_antibody_data)
#' check_data(example_antibody_data)
#' @export
check_data <- function(data,verbose=FALSE) {
    ## Check that all columns are present
    col.names <- c("individual", "sample_time", "biomarker_id","biomarker_group", "measurement", "repeat_number","birth")
    ## If there are any missing columns (NOTE: not checking if group or run are present)
    if (all(col.names %in% colnames(data)) != TRUE) {
        missing.cols <- col.names[which(col.names %in% colnames(data) == FALSE)] ## Find the missing column names
        if(verbose) message(paste(c("The following column(s) are missing from data: ", missing.cols), collapse = " "))

        ## Add missing variable names
        if(!("biomarker_group" %in% colnames(data))){
          if(verbose) message("Adding \"biomarker_group\" to data variables.")
          data$biomarker_group <- 1
        }
    
        if(!("repeat_number" %in% colnames(data))){
          if(verbose) message("Adding \"repeat_number\" to data variables.")
          data$repeat_number <- 1
        }
    }
    return(data)
}

#' Checks the entries of demographics data frame used for parameter stratification
#' @param demographics the data frame of demographics data. Must have columns: individual (integer ID of individual); birth (numeric time of birth); any other columns also in par_tab$stratification will be checked.
#' @param par_tab Optional: the parameter table controlling information such as bounds, initial values etc. If included, this will be used to check that the demographics data is consistent with the data.
#' @param verbose if TRUE, prints warning messages
#' @return the same data object with corrections if needed
#' @family check_inputs
#' @examples
#' data(example_antibody_data)
#' check_data(example_antibody_data)
#' @export
check_demographics <- function(demographics, par_tab=NULL, verbose=FALSE) {
  ## Check that all columns are present
  col.names <- c("individual", "birth")
  actual_colnames <- colnames(demographics)
  
  if (all(col.names %in% actual_colnames) != TRUE) {
    missing.cols <- col.names[which(col.names %in% actual_colnames == FALSE)] ## Find the missing column names
    if(verbose) message(paste(c("The following column(s) are missing from demographics: ", missing.cols), collapse = " "))
  }
  
  if(!is.null(par_tab)){
    strsplit1 <- function(x){
      if(!is.na(x)){
        return(strsplit(x,", "))
      } else {
        return(NA)
      }
    }
    
    ## Creates an estimated parameter entry for each 
    stratifications <- unique(unlist(sapply(par_tab[,"stratification"],function(x) strsplit1(x))))
    stratifications <- stratifications[!is.na(stratifications)]
    
    if (all(stratifications %in% actual_colnames) != TRUE) {
      missing.cols <- stratifications[which(stratifications %in% actual_colnames == FALSE)] ## Find the missing column names
      if(verbose) message(paste(c("The following column(s) are missing from demographics but requested in par_tab: ", stratifications), collapse = " "))
    }
  }  
  test_colnames <- actual_colnames[!(actual_colnames %in% c(col.names,"time","age"))]
  for(col in test_colnames){
    unique_levels <- unique(as.data.frame(demographics)[,col])
    if(min(unique_levels) != 0) warning(paste("Column ",col," in demographics should start at 0."))
    if(max(unique_levels) != length(unique_levels)-1) warning(paste("Column ",col," in demographics should be a sequence starting at 0."))
  }
  
  return(demographics)
}


#' Checks the attack_rates supplied in simulate_data
#' @param attack_rates a vector of attack_rates to be used in the simulation
#' @param possible_exposure_times vector of strain circulation times
#' @return nothing at the moment
#' @family check_inputs
#' @examples
#' attack_rates <- runif(40)
#' possible_exposure_times <- seq_len(40)
#' check_attack_rates(attack_rates, possible_exposure_times)
#' @export
check_attack_rates <- function(attack_rates, possible_exposure_times) {
    if(class(attack_rates) == "matrix") test_n <- ncol(attack_rates)
    else test_n <- length(attack_rates)
    if (test_n != length(possible_exposure_times)) stop("attack_rates is not the same length as possible_exposure_times")
    if (any(attack_rates < 0) || any(attack_rates > 1)) stop("attack_rates must be between 0 and 1")
}
#' Checks if the multivariate proposal is being used with the FOI proposal
#' @param version Which version of the posterior function to use? See \code{\link{create_posterior_func}}
#' @param mvr_pars Leave NULL to use univariate proposals. Otherwise, a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @return nothing at the moment
#' @family check_inputs
#' @export
check_proposals <- function(version, mvr_pars) {
    if (all(version == 1, !is.null(mvr_pars))) warning("The multivariate proposal can be inefficient for version 1.")
}

#' Check if the starting infection history table and titre data are consistent
#' @param antibody_data the data frame of titer data
#' @param possible_exposure_times the vector of times corresponding to entries in DOB
#' @param inf_hist the starting infection history matrix
#' @param verbose if TRUE, prints warning messages
#' @return nothing, prints a warning
#' @family check_inputs
#' @export
check_inf_hist <- function(antibody_data,possible_exposure_times, inf_hist,verbose=FALSE){
    DOBs <- create_age_mask(antibody_data %>% select(individual, birth) %>% distinct() %>% pull(birth),
                            possible_exposure_times)
    correct_dob <- rep(0,length(DOBs))
    for(i in seq_along(DOBs)){
        if(DOBs[i] > 1 & any(inf_hist[i,1:(DOBs[i]-1)] > 0)) correct_dob[i] <- 1
    }
    return(correct_dob)
}


