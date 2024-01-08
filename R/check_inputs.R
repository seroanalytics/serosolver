#' Check infection history matrix
#'
#' Checks that the infection history matrix is allowable given the birth dates and sampling times of the data
#' @param titre_dat the data frame of titre data
#' @param possible_exposure_times vector of times at which individuals could be exposed e.g., seq(1968,2015,by=1)
#' @param inf_hist the infection history matrix, with nrows = number indivs and ncols = length(possible_exposure_times)
#' @return a single boolean value, FALSE if the check passes, otherwise throws an error
#' @family check_inputs
#' data(example_titre_dat)
#' data(example_inf_hist)
#' times <- 1968:2015
#' check_inf_hist(example_titre_dat, times, example_inf_hist)
#' @export
check_inf_hist <- function(titre_dat, possible_exposure_times, inf_hist){
    DOBs <- get_DOBs(titre_dat)
    age_mask <- create_age_mask(DOBs[,2],possible_exposure_times)
    strain_mask <- create_strain_mask(titre_dat, possible_exposure_times)
    before_born <- logical(length(age_mask))
    after_sample <- logical(length(strain_mask))
    res <- logical(length(age_mask))
    for(i in seq_along(age_mask)){
        if(sum(inf_hist[i,] != 0)){
            first_inf <- min(which(inf_hist[i,] == 1))
            last_inf <- max(which(inf_hist[i,] == 1))
            
            before_born[i] <- first_inf < age_mask[i]
            after_sample[i] <- last_inf > strain_mask[i]
            res[i] <- before_born[i] | after_sample[i]
        } else {
            res[i] <- FALSE
        }
    }
    if (any(res)) {
        message(cat("Which infections before birth: ", which(before_born), "\n"))
        message(cat("Which infections after last sample: ", which(after_sample), "\n"))
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
#' @return the same par_tab object with corrections if needed
#' @family check_inputs
#' @examples
#' data(example_par_tab)
#' check_par_tab(example_par_tab, FALSE, version=1)
#' @export
check_par_tab <- function(par_tab, mcmc = FALSE, version = NULL) {
    ## Checks that should happen in simulate_data and run_MCMC
    essential_names <- c("names","values","fixed","steps","lower_bound","upper_bound","lower_start","upper_start","type")
    if (!all(essential_names %in% colnames(par_tab))) {
        message(paste(c("Some column names missing from par_tab: ", setdiff(essential_names,colnames(par_tab))),collapse=" "))
        if(!("steps" %in% colnames(par_tab))){
          message("Adding \"steps\" to par_tab variables.")
          par_tab$steps <- 0.1
        }
        if(!("type" %in% colnames(par_tab))){
          message("Adding \"obs_type\" to par_tab variables.")
          par_tab$obs_type <- 1
        }
    }
    pars <- par_tab$values
    names(pars) <- par_tab$names

    explicit_phi <- "phi" %in% names(pars)

    ## Additional checks that should happen in run_MCMC
    if (mcmc == TRUE) {
        phi_indices <- which(par_tab$names == "phi")
        no_phi <- length(phi_indices)
        explicit_phi <- (no_phi > 0)

        ## Check that all optional parameters are fixed, if not, fix them
        op_pars <- par_tab[which(par_tab$type == 0), ]
        if (all(op_pars$fixed == 1) == FALSE) stop("All optional parameters must be fixed")

        if (version == 1) {
            ## Check that the correct number of phis are present
            if(!explicit_phi){
                stop("prior version 1 is specified, but there are no entries for phi in par_tab. Note that there must be one phi entry per entry in possible_exposure_times")
            }
        }
        if (all(par_tab$fixed == 1)) {
            stop("all parameters are set to be fixed. The MCMC procedure will not work. Please set at least the entry for error to be fixed <- 0")
        }

        
        if (version %in% c(2, 3, 4)) {
            if (explicit_phi) stop(paste("phis are not required for versions 2, 3 or 4 but par_tab contains phi(s)")) ## Should we remove them?
        }
        ## Check bounds are equal to starting bounds
        if (any(par_tab$upper_start > par_tab$upper_bound) | any(par_tab$lower_start < par_tab$lower_bound)) {
            warning("Lower and upper bounds are not equal to the starting upper and lower bounds. If par_tab was used to create starting values, starting values may be out of bounds. ")
        }
    }
    ## Check that alpha and beta there for beta distribution
    ## If there, Pull out alpha and beta for beta binomial proposals
    if (!("alpha" %in% par_tab$names) | !("beta" %in% par_tab$names)) {
        stop("par_tab must have entries for `alpha` and `beta` for infection history prior")
    }
    par_tab
}

#' Checks the entries of data used in run_MCMC
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @return the same data object with corrections if needed
#' @family check_inputs
#' @examples
#' data(example_titre_dat)
#' check_data(example_titre_dat)
#' @export
check_data <- function(data) {
    ## Check that all columns are present
    col.names <- c("individual", "samples", "virus", "titre", "obs_type","group")
    ## If there are any missing columns (NOTE: not checking if group or run are present)
    if (all(col.names %in% colnames(data)) != TRUE) {
        missing.cols <- col.names[which(col.names %in% colnames(data) == FALSE)] ## Find the missing column names
        message(paste(c("The following column(s) are missing from data: ", missing.cols), collapse = " "))

        if(!("type" %in% colnames(data))){
          message("Adding \"obs_type\" to data variables.")
          data$obs_type <- 1
        }
    }
    return(data)
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
    if (length(attack_rates) != length(possible_exposure_times)) stop("attack_rates is not the same length as possible_exposure_times")
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
#' @param titre_dat the data frame of titer data
#' @param possible_exposure_times the vector of times corresponding to entries in DOB
#' @param inf_hist the starting infection history matrix
#' @return nothing, prints a warning
#' @family check_inputs
#' @export
check_inf_hist <- function(titre_dat,possible_exposure_times, inf_hist){
    DOBs <- create_age_mask(titre_dat %>% select(individual, DOB) %>% distinct() %>% pull(DOB),
                            possible_exposure_times)
    correct_dob <- rep(0,length(DOBs))
    for(i in seq_along(DOBs)){
        if(DOBs[i] > 1 & any(inf_hist[i,1:(DOBs[i]-1)] > 0)) correct_dob[i] <- 1
    }
    return(correct_dob)
}
