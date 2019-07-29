#' Check par_tab for simulate_data
#'
#' Checks the entries of par_tab used in simulate_data
#' @param par_tab the parameter table controlling information such as bounds, initial values etc
#' @param mcmc logical, if TRUE then checks are performed for the MCMC algorithm. Use FALSE when simulating data
#' @param version which version of the posterior function is being used? See \code{\link{create_posterior_func}}
#' @return nothing, only an error message if necessary
#' @family check_inputs
#' @examples
#' data(example_par_tab)
#' check_par_tab(example_par_tab, FALSE, version=1)
#' @export
check_par_tab <- function(par_tab, mcmc = FALSE, version = NULL) {
    ## Checks that should happen in simulate_data and run_MCMC
    essential_names <- c("names","values","fixed","steps","lower_bound","upper_bound","lower_start","upper_start","type")
    if (!all(essential_names %in% colnames(par_tab))) {
        stop(paste(c("Some column names missing from par_tab: ", setdiff(essential_names,colnames(par_tab))),collapse=" "))
    }
    pars <- par_tab$values
    names(pars) <- par_tab$names

    ## Checks for wane_type
    ## Extract wane_type
    wane_type <- pars["wane_type"]
    if (is.na(wane_type)) stop("wane_type is missing from par_tab (specify 0 for linear decrease or 1 for piecewise linear) ") ## If user has not entered wane_type in par_tab
                                        # Check that additional parameters are present
    if (wane_type == 1 & (!("kappa" %in% names(pars)) | !("t_change" %in% names(pars)))) stop("Parameters needed for wane_type=1 (piecewise linear) are missing")
    explicit_phi <- "phi" %in% names(pars)

    ## Additional checks that should happen in run_MCMC
    if (mcmc == TRUE) {
        phi_indices <- which(par_tab$type == 2)
        no_phi <- length(phi_indices)
        explicit_phi <- (no_phi > 0)

        ## Check that all optional parameters are fixed, if not, fix them
        op_pars <- par_tab[which(par_tab$type == 0), ]
        if (all(op_pars$fixed == 1) == FALSE) stop("All optional parameters must be fixed")

        if (version == 1) {
            ## Check that the correct number of phis are present
            if(!explicit_phi){
                stop("prior version 1 is specified, but there are no entries for phi in par_tab. Note that there must be one phi entry per entry in strain_isolation_times")
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

}

#' Checks the entries of data used in run_MCMC
#' @param data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @return nothing at the moment
#' @family check_inputs
#' @examples
#' data(example_titre_dat)
#' check_data(example_titre_dat)
#' @export
check_data <- function(data) {
    ## Check that all columns are present
    col.names <- c("individual", "samples", "virus", "titre", "group")
    ## If there are any missing columns (NOTE: not checking if group or run are present)
    if (all(col.names %in% colnames(data)) != TRUE) {
        missing.cols <- col.names[which(col.names %in% colnames(data) == FALSE)] ## Find the missing column names
        stop(paste(c("The following column(s) are missing from data: ", missing.cols), collapse = " "))
    }
}

#' Checks the attack_rates supplied in simulate_data
#' @param attack_rates a vector of attack_rates to be used in the simulation
#' @param strain_isolation_times vector of strain circulation times
#' @return nothing at the moment
#' @family check_inputs
#' @examples
#' attack_rates <- runif(40)
#' strain_isolation_times <- seq_len(40)
#' check_attack_rates(attack_rates, strain_isolation_times)
#' @export
check_attack_rates <- function(attack_rates, strain_isolation_times) {
    if (length(attack_rates) != length(strain_isolation_times)) stop("attack_rates is not the same length as strain_isolation_times")
    if (any(attack_rates < 0) || any(attack_rates > 1)) stop("attack_rates must be between 0 and 1")
}
#' Checks if the multivarite proposal is being used with the FOI proposal
#' @param version Which version of the posterior function to use? See \code{\link{create_post_func}}
#' @param mvr_pars Leave NULL to use univariate proposals. Otherwise, a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @return nothing at the moment
#' @family check_inputs
#' @export
check_proposals <- function(version, mvr_pars) {
    if (all(version == 1, !is.null(mvr_pars))) warning("The multivariate proposal can be inefficient for version 1.")
}
