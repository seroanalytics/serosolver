#' Likelihood function given discrete data (normal)
#'
#' Calculates the likelihood of observing a set of discrete HI titres given a corresponding set of predicted titres
#' @param expected vector of expected HI titres
#' @param data vector of observed discrete HI titres
#' @param expected_indices the indices of the measurement_shifts vector that each predicted titre needs adding to it
#' @param measurement_shifts the vector of measurement shifts for each cluster to add to the predicted titres
#' @return a vector with the likelihood of making each observation given the predictions
#' @export
r_likelihood <- function(expected, data, theta, expected_indices = NULL, measurement_shifts = NULL) {
  if (!is.null(expected_indices) & !is.null(measurement_shifts)) {
    expected <- expected + measurement_shifts[expected_indices]
  }

  ## Vectorise, calculate boundaries separately
  liks <- numeric(length(expected))
  large_i <- data >= theta["MAX_TITRE"]
  small_i <- data < 1
  rest_i <- data >= 1 & data < theta["MAX_TITRE"]

  liks[large_i] <- pnorm(theta["MAX_TITRE"], expected[large_i], theta["error"], lower.tail = FALSE, log.p = TRUE)
  liks[small_i] <- pnorm(1, expected[small_i], theta["error"], lower.tail = TRUE, log.p = TRUE)
  liks[rest_i] <- log(pnorm(data[rest_i] + 1, expected[rest_i], theta["error"], lower.tail = TRUE, log.p = FALSE) -
    pnorm(data[rest_i], expected[rest_i], theta["error"], lower.tail = TRUE, log.p = FALSE))
  return(liks)
}


#' Likelihood function given continuous data (normal)
#'
#' Calculates the likelihood of observing a set of continuous and bounded titres given a corresponding set of predicted titres
#' @param expected vector of expected titres
#' @param data vector of observed continuous titres
#' @param expected_indices the indices of the measurement_shifts vector that each predicted titre needs adding to it
#' @param measurement_shifts the vector of measurement shifts for each cluster to add to the predicted titres
#' @return a vector with the likelihood of making each observation given the predictions
#' @export
r_likelihood_continuous <- function(expected, data, theta, expected_indices = NULL, measurement_shifts = NULL) {
  if (!is.null(expected_indices) & !is.null(measurement_shifts)) {
    expected <- expected + measurement_shifts[expected_indices]
  }
  
  ## Vectorise, calculate boundaries separately
  liks <- numeric(length(expected))
  large_i <- data >= theta["MAX_TITRE"]
  small_i <- data <= theta["MIN_TITRE"]
  rest_i <- data > theta["MIN_TITRE"] & data < theta["MAX_TITRE"]
  
  liks[large_i] <- pnorm(theta["MAX_TITRE"], expected[large_i], theta["error"], lower.tail = FALSE, log.p = TRUE)
  liks[small_i] <- pnorm(theta["MIN_TITRE"], expected[small_i], theta["error"], lower.tail = TRUE, log.p = TRUE)
  liks[rest_i] <- dnorm(data[rest_i], expected[rest_i], theta["error"], log = TRUE)
  return(liks)
}



#' Calculate FOI log probability
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @param phis a vector of FOIs
#' @param infection_history the matrix of infection histories
#' @param age_mask the age mask vector as returned by \code{\link{create_age_mask}}
#' @return a single log probability
#' @family priors
#' @export
calc_phi_probs <- function(phis, infection_history, age_mask, strain_mask) {
  lik <- 0
  for (i in 1:ncol(infection_history)) {
    use_indivs <- intersect(which(age_mask <= i), which(strain_mask >= i))
    lik <- lik + sum(log((phis[i]^infection_history[use_indivs, i] *
      (1 - phis[i])^(1 - infection_history[use_indivs, i]))))
  }
  lik
}

#' Calculate FOI log probability vector
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @inheritParams calc_phi_probs
#' @return a vector of log probabilities for each individual
#' @family priors
#' @export
calc_phi_probs_indiv <- function(phis, infection_history, age_mask, strain_mask) {
  lik <- numeric(nrow(infection_history))
  for (i in 1:ncol(infection_history)) {
    to_add <- log(((phis[i]^infection_history[, i]) *
      (1 - phis[i])^(1 - infection_history[, i]))) *
      as.numeric(age_mask <= i) *
      as.numeric(strain_mask >= i)
    lik <- lik + to_add
  }
  lik
}

calc_phi_loc_probs_indiv <- function(phis, group_probs, infection_history, age_mask, strain_mask, group_indices) {
  lik <- numeric(nrow(infection_history))
  max_group_p <- max(group_probs)
  group_probs <- group_probs / max_group_p
  for (i in 1:ncol(infection_history)) {
    tmp <- (infection_history[, i] * log(phis[i] * group_probs[group_indices]) +
      (1 - infection_history[, i]) * log(1 - phis[i] * group_probs[group_indices])) *
      as.numeric(age_mask <= i) * as.numeric(strain_mask >= i)
    lik <- lik + tmp
  }
  lik
}


#' Calculate FOI from spline - INACTIVE
#'
#' Version of FOI prior that calculates a seasonal spline such that the FOI in a given time point (month) comes from this spline term
#' @inheritParams calc_phi_probs
#' @param foi vector of force of infections per year
#' @param knots vector of knots for the spline
#' @param theta vector of theta parameters for spline
#' @return a vector of log probabilities for each individual
#' @family priors
calc_phi_probs_spline <- function(foi, knots, theta, infection_history, age_mask) {
  phis <- generate_phis(foi, knots, theta, length(foi), 12)
  lik <- numeric(nrow(infection_history))
  for (i in 1:ncol(infection_history)) {
    lik <- lik + infection_history[, i] * log(phis[i]) + (1 - infection_history[, i]) * log(phis[i]) + log(as.numeric(age_mask <= i)) + log(as.numeric(strain_mask >= i))
  }
  lik
}

#' Generate FOI phis - INACTIVE
#'
#' Calculates the seasonal spline for \code{\link{calc_phi_probs_spline}}
#' @inheritParams calc_phi_probs_spline
#' @param n_years number of years to calculate
#' @param buckets number of buckets per year (12 for monthly, 1 for annual)
#' @param degree degree of the spline
#' @return a vector of FOIs for each time point
generate_phis <- function(foi, knots, theta, n_years, buckets, degree = 2) {
  x <- seq(0, buckets - 1, by = 1) / buckets
  n_knots <- length(knots) + degree + 1
  all_dat <- NULL
  index <- 1
  tmp <- gen_spline_y(x, knots, degree, theta)
  tmp <- tmp / sum(tmp)
  all_dat <- numeric(length(tmp) * length(foi))
  for (i in 1:n_years) {
    all_dat[index:(index + length(tmp) - 1)] <- tmp * foi[i]
    index <- index + length(tmp)
  }
  all_dat
}

#' Generates a spline for \code{\link{generate_phis}} - INACTIVE
gen_spline_y <- function(x, knots, degree, theta, intercept = TRUE) {
  basis <- bs(
    x = x, knots = knots, degree = degree,
    Boundary.knots = c(0, 1), intercept = intercept
  )

  y.spline <- basis %*% theta
  return(as.vector(y.spline))
}


#' Prior on measurement shifts
#'
#' Assumes measurement shifts are drawn from a normal distribution (random effects) with given standard deviation and mean. Code is commented out to assume log normally distributied
#' @param rhos vector of measurement shifts
#' @param pars vector of parameters, including rho_mean and rho_sd for the normal distribution
#' @return a single prior probability
#' @family priors
#' @export
prob_shifts <- function(rhos, pars) {
  rho_mean <- pars["rho_mean"]
  rho_sd <- pars["rho_sd"]
  ## l_mean <- log(mu_mean) - (mu_sd^2)/2
  return(sum(dnorm(rhos, rho_mean, rho_sd, log = TRUE)))
  # return(sum(dlnorm(mus, l_mean, mu_sd, log=TRUE)))
}

#' Measurement shift creation
#'
#' Creates a function to calculate the prior probability of a set of measurement shifts. Assumes normally distributed random effects
#' @param par_tab the parameter table as in \code{\link{create_posterior_func}}
#' @return a function pointer to solve the measurement shifts prior
#' @family priors
#' @export
create_prob_shifts <- function(par_tab) {
  par_names <- par_tab$names
  rho_indices <- which(par_tab$type == 3)

  f <- function(pars) {
    names(pars) <- par_names
    rhos <- pars[rho_indices]
    rho_mean <- pars["rho_mean"]
    rho_sd <- pars["rho_sd"]
    return(sum(dnorm(rhos, rho_mean, rho_sd, log = TRUE)))
  }
  f
}

#' Create strain specific bias prior
#'
#' Creates a function that can be used to easily solve the prior on the strain-specific boosting parameters (normal random effects)
#' @inheritParams create_prob_shifts
#' @return a function pointer
#' @family priors
#' @export
create_prior_mu <- function(par_tab) {
  ## Extract parameter type indices from par_tab, to split up
  ## similar parameters in model solving functions
  option_indices <- which(par_tab$type == 0)
  theta_indices <- which(par_tab$type %in% c(0, 1))
  phi_indices <- which(par_tab$type == 2)
  measurement_indices_par_tab <- which(par_tab$type == 3)
  weights_indices <- which(par_tab$type == 4) ## For functional form version
  knot_indices <- which(par_tab$type == 5)
  mu_indices_par_tab <- which(par_tab$type == 6)

  par_names_theta <- par_tab[theta_indices, "names"]

  ## Expect to only take the vector of parameters
  f <- function(pars) {
    mus <- pars[mu_indices_par_tab]
    pars <- pars[theta_indices]
    names(pars) <- par_names_theta
    return(prob_mus(mus, pars))
  }
}

#' Prior probability of strain specific boosting
#'
#' Function find the prior probability of a set of boosting parameters assuming that boosting is drawn from a normal distribution (random effects)
#' @param mus the vector of boosting parameters
#' @param pars the vector of other model parameters including mu_mean and mu_sd
#' @return a single log prior probability
#' @family priors
#' @export
prob_mus <- function(mus, pars) {
  mu_mean <- pars["mu_mean"]
  mu_sd <- pars["mu_sd"]
  return(sum(dnorm(mus, mu_mean, mu_sd, log = TRUE)))
  ## 
  ## location <- log(mu_mean^2 / sqrt(mu_sd^2 + mu_mean^2))
  ## shape <- sqrt(log(1 + (mu_sd^2/mu_mean^2)))
  ## l_mean <- log(mu_mean) - (mu_sd^2)/2
  ## p <- sum(dnorm(log(mus),mu_mean,mu_sd,log=TRUE))
  ## p_mean <- 0.6
  ## p_sd <- 0.5
  ## p_mu <- log(p_mean/sqrt(1 + (p_sd/p_mean)^2))
  ## p_sigma <- sqrt(log(1 + (p_sd/p_mean)^2))
  ## p_lik <- log(p_sigma*2.506628) - 0.5*((mu_mean - p_mu)/p_sigma)^2
  ## return(p+p_lik)
  ## return(sum(log(dtruncnorm(mus, a=0,mean=mu_mean, sd=mu_sd))))
  ## return(sum(dlnorm(mus, location, shape, log=TRUE)))
  ## mean_log_y <- mean(log(mus))
  ## sd_log_y <- sd(log(mus))
  ## sigmaOfLogY <- dunif(mu_sd, 0.001*sd_log_y,1000*sd_log_y)
  ## muOfLogY <- dnorm(mu_mean, mean_log_y, 1/(10*sd_log_y)^2)
  ## return(sum(dlnorm(mus, mu_mean, 1/mu_sd^2, log=TRUE)) + sigmaOfLogY + muOfLogY)
  ## return(sum(dlnorm(mus, mu_mean, mu_sd, log = TRUE)))
}


#' Posterior function pointer
#'
#' Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to pass the data set back and forth. No example is provided for function_type=2, as this should only be called within \code{\link{run_MCMC}}
#' @param par_tab the parameter table controlling information such as bounds, initial values etc. See \code{\link{example_par_tab}}
#' @param titre_dat the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); obs_type (integer of the observation group type, using a unique value for each distinctive type of observation underpinned by the same generative model); titre (integer of titre value against the given virus at that sampling time). See \code{\link{example_titre_dat}}
#' @param antigenic_map (optional) a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. See \code{\link{example_antigenic_map}}
#' @param strain_isolation_times (optional) if no antigenic map is specified, this argument gives the vector of times at which individuals can be infected
#' @param version which infection history assumption version to use? See \code{\link{describe_priors}} for options. Can be 1, 2, 3 or 4
#' @param solve_likelihood usually set to TRUE. If FALSE, does not solve the likelihood and instead just samples/solves based on the model prior
#' @param age_mask see \code{\link{create_age_mask}} - a vector with one entry for each individual specifying the first epoch of circulation in which an individual could have been exposed
#' @param measurement_indices_by_time if not NULL, then use these indices to specify which measurement bias parameter index corresponds to which time
#' @param mu_indices if not NULL, then use these indices to specify which boosting parameter index corresponds to which time
#' @param n_alive if not NULL, uses this as the number alive in a given year rather than calculating from the ages. This is needed if the number of alive individuals is known, but individual birth dates are not
#' @param function_type integer specifying which version of this function to use. Specify 1 to give a posterior solving function; 2 to give the gibbs sampler for infection history proposals; otherwise just solves the titre model and returns predicted titres. NOTE that this is not the same as the attack rate prior argument, \code{version}!
#' @param titre_before_infection TRUE/FALSE value. If TRUE, solves titre predictions, but gives the predicted titre at a given time point BEFORE any infection during that time occurs.
#' @param data_type integer, currently accepting 1 or 2. Set to 1 for discretized, bounded data, or 2 for continuous, bounded data. Note that with 2, MIN_TITRE must be set.
#' @param ... other arguments to pass to the posterior solving function
#' @return a single function pointer that takes only pars and infection_histories as unnamed arguments. This function goes on to return a vector of posterior values for each individual
#' @examples
#' \dontrun{
#' data(example_par_tab)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' data(example_inf_hist)
#'
#' ## Simple model solving code. Output matches entries of example_titre_dat
#' model_func <- create_posterior_func(example_par_tab, example_titre_dat, example_antigenic_map, function_type = 3)
#' y <- model_func(example_par_tab$values, example_inf_hist)
#'
#' ## Solve likelihood
#' par_tab <- example_par_tab[example_par_tab$names != "phi",]
#' likelihood_func <- create_posterior_func(par_tab, example_titre_dat, example_antigenic_map, function_type = 1, version = 2)
#' liks <- likelihood_func(par_tab$values, example_inf_hist)
#' }
#' @export
create_posterior_func <- function(par_tab,
                                  titre_dat,
                                  antigenic_map=NULL,
                                  strain_isolation_times=NULL,
                                  version = 1,
                                  solve_likelihood = TRUE,
                                  age_mask = NULL,
                                  measurement_indices_by_time = NULL,
                                  mu_indices = NULL,
                                  n_alive = NULL,
                                  function_type = 1,
                                  titre_before_infection=FALSE,
                                  data_type=1,
                                  obs_types_weights =1,
                                  ...) {
    check_par_tab(par_tab, TRUE, version)
    if (!("group" %in% colnames(titre_dat))) {
        titre_dat$group <- 1
    }
   
    ## Add a dummy observation type variable if not provided
    if (!("obs_type" %in% colnames(titre_dat))) {
        message(cat("Note: no obs_type detection in titre_dat. Assuming all obs_type as 1."))
        titre_dat$obs_type <- 1
    }
    
    if (!("obs_type" %in% colnames(par_tab))) {
        message(cat("Note: no obs_type detection in par_tab Assuming all obs_type as 1."))
        par_tab$obs_type <- 1
    }
    
    check_data(titre_dat)
    
    titre_dat <- titre_dat %>% arrange(individual, obs_type, samples, virus, run)
    
    ## Get unique observation types
    unique_obs_types <- unique(titre_dat$obs_type)
    n_obs_types <- length(unique_obs_types)
    
    ## Likelihood versions for different obs types
    if(length(data_type) ==1 & n_obs_types > 1){
        data_type <- rep(data_type, n_obs_types)
    }
    
    if(length(obs_types_weights) ==1 & n_obs_types > 1){
        obs_types_weights <- rep(1, n_obs_types)
    }
    
    #########################################################
    ## SETUP ANTIGENIC MAP
    #########################################################
    ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
    if (!is.null(antigenic_map)) {
        strain_isolation_times_tmp <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
        if(!is.null(strain_isolation_times) & !identical(strain_isolation_times, strain_isolation_times_tmp)){
            message(cat("Warning: provided strain_isolation_times argument does not match entries in the antigenic map. Please make sure that there is an entry in the antigenic map for each possible circulation time. Using the antigenic map times."))
        }
      strain_isolation_times <- strain_isolation_times_tmp
      
      ## If no observation types assumed, set all to 1.
      if (!("obs_type" %in% colnames(antigenic_map))) {
          message(cat("Note: no obs_type detection in antigenic_map. Aligning antigenic map with par_tab."))
          antigenic_map_tmp <- replicate(n_obs_types,antigenic_map,simplify=FALSE)
          for(obs_type in unique_obs_types){
              antigenic_map_tmp[[obs_type]]$obs_type <- obs_type
          }
          antigenic_map <- do.call(rbind,antigenic_map_tmp)
      }
      
    } else {
        ## Create a dummy map with entries for each observation type
      antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,
                                  "inf_times"=rep(strain_isolation_times, n_obs_types), 
                                  "obs_type"=rep(unique_obs_types,each=length(strain_isolation_times)))
    }
    #########################################################
    ## SETUP DATA
    #########################################################
    ## Separate out initial readings and repeat readings - we only
    ## want to solve the model once for each unique indiv/sample/virus year tested
    titre_dat_unique <- titre_dat[titre_dat$run == 1, ]
    ## Observations from repeats
    titre_dat_repeats <- titre_dat[titre_dat$run != 1, ]
    ## Find which entry in titre_dat_unique each titre_dat_repeats entry should correspond to
    tmp <- row.match(
        titre_dat_repeats[, c("individual", "samples", "obs_type", "virus")],
        titre_dat_unique[, c("individual", "samples", "obs_type", "virus")]
    )
    titre_dat_repeats$index <- tmp

    ## Which entries in the overall titre_dat matrix does each entry in titre_dat_unique correspond to?
    overall_indices <- row.match(
        titre_dat[, c("individual", "samples", "obs_type","virus")],
        titre_dat_unique[, c("individual", "samples", "obs_type","virus")]
    )

    ## Setup data vectors and extract
    setup_dat <- setup_titredat_for_posterior_func(
        titre_dat_unique, antigenic_map, 
        strain_isolation_times,
        age_mask, n_alive
    )
    ## Vector of observation types matching the unique samples
    obs_types <- setup_dat$obs_types
    
    ## Number of unique groups
    n_groups <- length(unique(titre_dat$group))
    group_id_vec <- setup_dat$group_id_vec
    
    ## List of melted antigenic maps, one entry for each observation type
    antigenic_map_melted <- setup_dat$antigenic_map_melted
    antigenic_distances <- antigenic_map_melted[[1]]
    
    strain_isolation_times <- setup_dat$strain_isolation_times
    infection_strain_indices <- setup_dat$infection_strain_indices
    
    ## Sample collection times, entry for each unique individual, observation type and sample
    sample_times <- setup_dat$sample_times
    ## Indices related to entries in sample_data
    sample_data_start <- setup_dat$sample_data_start
    
    ## Indices related to entries in titre_dat
    nrows_per_sample <- setup_dat$nrows_per_sample
    titre_data_start <- setup_dat$titre_data_start
    
    ## Indices related to entries in type_data
    type_data_start <- setup_dat$type_data_start
    obs_types <- setup_dat$obs_types
    
    ## Indices related to entries in the antigenic map
    measured_strain_indices <- setup_dat$measured_strain_indices
    
    n_alive <- setup_dat$n_alive
    age_mask <- setup_dat$age_mask
    strain_mask <- setup_dat$strain_mask
    n_indiv <- setup_dat$n_indiv
    DOBs <- setup_dat$DOBs

    ## Which entries of the unique titre data correspond to each individual? 
    ## Used to summarize into per-individual likelihoods later
    nrows_per_individual_in_data <- lapply(unique_obs_types, function(y) 
        plyr::ddply(titre_dat_unique[titre_dat_unique$obs_type == y,], .(individual), "nrow")$nrow
        )
    cum_nrows_per_individual_in_data <- lapply(nrows_per_individual_in_data, function(x) cumsum(c(0, x)))
    
    ## Some additional setup for the repeat data
    ## Used to summarize into per-individual likelihoods later
    nrows_per_individual_in_data_repeats <- lapply(unique_obs_types, function(y) 
        plyr::ddply(titre_dat_repeats[titre_dat_repeats$obs_type == y,], .(individual), "nrow")$nrow
    )
    cum_nrows_per_individual_in_data_repeats <- lapply(nrows_per_individual_in_data_repeats, function(x) cumsum(c(0, x)))
    
    obs_type_indices <- lapply(unique_obs_types, function(x) which(titre_dat_unique$obs_type == x))
    obs_type_indices_repeats <- lapply(unique_obs_types, function(x) which(titre_dat_repeats$obs_type == x))
    
    ## Pull out unique and repeat titres for solving likelihood later
    titres_unique <- titre_dat_unique$titre
    titres_repeats <- titre_dat_repeats$titre
    repeat_indices <- titre_dat_repeats$index
    repeat_indices_cpp <- repeat_indices - 1
    
    #########################################################
    ## PARAMETER TABLE
    #########################################################
    ## Extract parameter type indices from par_tab, to split up
    ## similar parameters in model solving functions
    
    ## In general we are just going to use the indices for a single observation type
    par_tab_unique <- par_tab[!is.na(par_tab$obs_type) & par_tab$obs_type == min(par_tab$obs_type),]
    theta_indices_unique <- which(par_tab_unique$type %in% c(0, 1))
    
    ## Each obs_type must have the same vector of parameters in the same order
    par_names_theta <- par_tab_unique[theta_indices_unique, "names"]
    theta_indices_unique <- theta_indices_unique - 1
    names(theta_indices_unique) <- par_names_theta
    n_pars <- length(theta_indices_unique)

    ## These will be different for each obs_type
    option_indices <- which(par_tab$type == 0)
    theta_indices <- which(par_tab$type %in% c(0, 1))
    measurement_indices_par_tab <- which(par_tab$type == 3)
    mu_indices_par_tab <- which(par_tab$type == 6)
    
    par_names_theta_all <- par_tab[theta_indices,"names"]
    
    
    ## Sort out any assumptions for measurement bias
    use_measurement_bias <- (length(measurement_indices_par_tab) > 0) & !is.null(measurement_indices_by_time)
    
    titre_shifts <- c(0)
    expected_indices <- NULL
    measurement_bias <- NULL
    use_strain_dependent <- (length(mu_indices) > 0) & !is.null(mu_indices)
    additional_arguments <- NULL

    repeat_data_exist <- nrow(titre_dat_repeats) > 0

    if (use_measurement_bias) {
        message(cat("Using measurement bias\n"))
        expected_indices <- measurement_indices_by_time[match(titre_dat_unique$virus, strain_isolation_times)]
    } else {
        expected_indices <- c(-1)
    }

    if (use_strain_dependent) {
        boosting_vec_indices <- mu_indices - 1
        mus <- rep(2, length(strain_isolation_times))
    } else {
        boosting_vec_indices <- mus <- c(-1)
    }

    if (!repeat_data_exist) {
        repeat_indices_cpp <- c(-1)
    }

    ## These will be the same for each obs_type, as currently only one exposure type
    phi_indices <- which(par_tab$type == 2)
    weights_indices <- which(par_tab$type == 4) ## For functional form version of FOI
    knot_indices <- which(par_tab$type == 5) ## For functional form version of FOI
    
    ## Find which options are being used in advance for speed
    explicit_phi <- (length(phi_indices) > 0) ## Explicit prob of infection term
    spline_phi <- (length(knot_indices) > 0)
    
    ## Allow different likelihood function for each observation type
    ## Set data type for likelihood function
    ## Potentially different likelihood for each observation type
    likelihood_func_use <- list()
    for(obs_type in unique_obs_types){
        if(data_type[obs_type] == 1){
          likelihood_func_use[[obs_type]] <- likelihood_func_fast
          message(cat("Setting to discretized, bounded observations\n"))
          
        } else if(data_type[obs_type] == 2){
          message(cat("Setting to continuous, bounded observations\n"))
          likelihood_func_use[[obs_type]] <- likelihood_func_fast_continuous
        } else {
          message(cat("Assuming discretized, bounded observations\n"))
          likelihood_func_use[[obs_type]] <- likelihood_func_fast
        }
    }
    
    if (function_type == 1) {
        message(cat("Creating posterior solving function...\n"))
        f <- function(pars, infection_history_mat) {
            theta <- pars[theta_indices]
            names(theta) <- par_names_theta_all
            
            ## Pass strain-dependent boosting down
            if (use_strain_dependent) {
                mus <- pars[mu_indices_par_tab]
            }
            
            antigenic_map_long <- matrix(nrow=length(strain_isolation_times)^2, ncol=n_obs_types)
            antigenic_map_short <- matrix(nrow=length(strain_isolation_times)^2, ncol=n_obs_types)
            
            sigma1s <- theta[which(par_names_theta_all=="sigma1")]
            sigma2s <- theta[which(par_names_theta_all=="sigma2")]
            
            for(obs_type in unique_obs_types){
                antigenic_map_long[,obs_type] <- create_cross_reactivity_vector(antigenic_map_melted[[obs_type]], sigma1s[obs_type])
                antigenic_map_short[,obs_type] <- create_cross_reactivity_vector(antigenic_map_melted[[obs_type]], sigma2s[obs_type])
            }
            
            y_new <- titre_data_fast(
                theta, theta_indices_unique, unique_obs_types,
                infection_history_mat, strain_isolation_times, infection_strain_indices,
                sample_times, type_data_start,obs_types,
                sample_data_start, titre_data_start,
                nrows_per_sample, measured_strain_indices, 
                antigenic_map_long,
                antigenic_map_short,
                antigenic_distances,
                mus, boosting_vec_indices,
                titre_before_infection
            )
            
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_par_tab]
                titre_shifts <- measurement_bias[expected_indices]
                y_new <- y_new + titre_shifts
            }
            
            ## Transmission prob is the part of the likelihood function corresponding to each individual
            transmission_prob <- rep(0, n_indiv)
            if (explicit_phi) {
                phis <- pars[phi_indices]
                if (spline_phi) {
                    ## If using spline term for FOI, add here
                    weights <- pars[weights_indices]
                    knots <- pars[knot_indices]
                    liks <- liks + calc_phi_probs_spline(
                                       phis, knots, weights,
                                       infection_history_mat, age_mask
                                   )
                } else {
                    ## Or the baseline transmission likelihood contribution
                    transmission_prob <- calc_phi_probs_indiv(
                        phis, infection_history_mat,
                        age_mask, strain_mask
                    )
                }
            }
            if (solve_likelihood) {
                ## Calculate likelihood for unique titres and repeat data
                ## Sum these for each individual
                liks <- numeric(n_indivs)
                for(obs_type in unique_obs_types){
                    ## Need theta for each observation type
                    liks_tmp <- likelihood_func_use[[obs_type]](
                                                    theta[(theta_indices_unique+1) + n_pars*(obs_type-1)], 
                                                    titres_unique[obs_type_indices[[obs_type]]], 
                                                    y_new[obs_type_indices[[obs_type]]])
                    
                    liks <- liks + obs_types_weights[obs_type]*sum_buckets(liks_tmp, nrows_per_individual_in_data[[obs_type]])
                    if (repeat_data_exist) {
                        ## Need theta for each observation type
                        
                        liks_repeats <- likelihood_func_use[[obs_type]](
                            theta[(theta_indices_unique+1) + n_pars*(obs_type-1)], 
                            titres_repeats[obs_type_indices_repeats[[obs_type]]], 
                            y_new[repeat_indices][obs_type_indices_repeats[[obs_type]]])
                        
                        liks <- liks + obs_types_weights[obs_type]*sum_buckets(liks_repeats, nrows_per_individual_in_data_repeats[[obs_type]])
                    }
                }
            } else {
                liks <- rep(-100000, n_indiv)
            }
            return(list(liks, transmission_prob))
        }
    } else if (function_type == 2) {
        
        message(cat("Creating infection history proposal function\n"))
        if (version == 4) {
            n_alive_total <- rowSums(n_alive)
        } else {
            n_alive_total <- c(-1, -1)
        }
        alpha <- par_tab[par_tab$names == "alpha","values"]
        beta <- par_tab[par_tab$names == "beta","values"]
        n_infected_group <- c(0, 0)
        ## Generate prior lookup table
        lookup_tab <- create_prior_lookup_groups(titre_dat, strain_isolation_times, alpha, beta, n_alive)
        ## Use the original gibbs proposal function if no titre immunity
        f <- function(pars, infection_history_mat,
                      probs, sampled_indivs,
                      alpha, beta,
                      n_infs, swap_propn,
                      swap_dist,
                      proposal_iter,
                      accepted_iter,
                      proposal_swap,
                      accepted_swap,
                      overall_swap_proposals,
                      overall_add_proposals,
                      proposal_ratios,
                      temp=1,
                      propose_from_prior=TRUE) {
            theta <- pars[theta_indices]
            names(theta) <- par_names_theta

            ## Pass strain-dependent boosting down
            if (use_strain_dependent) {
                mus <- pars[mu_indices_par_tab]
            }
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_par_tab]
                titre_shifts <- measurement_bias[expected_indices]
            }
            ## Work out short and long term boosting cross reactivity - C++ function
            antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])

            n_infections <- sum_infections_by_group(infection_history_mat, group_id_vec, n_groups)
            if (version == 4) n_infected_group <- rowSums(n_infections)
            ## Now pass to the C++ function
            res <- inf_hist_prop_prior_v2_and_v4(
                theta,
                infection_history_mat,
                probs,
                sampled_indivs,
                n_infs,
                age_mask,
                strain_mask,
                n_alive,
                n_infections,
                n_infected_group,
                lookup_tab,
                swap_propn,
                swap_dist,
                propose_from_prior,
                alpha,
                beta,
                strain_isolation_times,
                infection_strain_indices,
                sample_times,
                rows_per_indiv_in_samples,
                cum_nrows_per_individual_in_data,
                cum_nrows_per_individual_in_data_repeats,
                nrows_per_sample,
                group_id_vec,
                measured_strain_indices,
                antigenic_map_long,
                antigenic_map_short,
                antigenic_distances,
                titres_unique,
                titres_repeats,
                repeat_indices_cpp,
                titre_shifts,
                proposal_iter = proposal_iter,
                accepted_iter = accepted_iter,
                proposal_swap = proposal_swap,
                accepted_swap = accepted_swap,
                overall_swap_proposals,
                overall_add_proposals,
                proposal_ratios,
                mus,
                boosting_vec_indices,
                n_alive_total,
                temp,
                solve_likelihood,
                data_type
            )
            return(res)
        }
    } else {
        message(cat("Creating model solving function...\n"))
        ## Final version is just the model solving function
        f <- function(pars, infection_history_mat) {
            theta <- pars[theta_indices]
            names(theta) <- par_names_theta_all

            ## Pass strain-dependent boosting down
            if (use_strain_dependent) {
                mus <- pars[mu_indices_par_tab]
            }
            
            antigenic_map_long <- matrix(nrow=length(strain_isolation_times)^2, ncol=n_obs_types)
            antigenic_map_short <- matrix(nrow=length(strain_isolation_times)^2, ncol=n_obs_types)
            
            sigma1s <- theta[which(par_names_theta_all=="sigma1")]
            sigma2s <- theta[which(par_names_theta_all=="sigma2")]
            
            for(obs_type in unique_obs_types){
                antigenic_map_long[,obs_type] <- create_cross_reactivity_vector(antigenic_map_melted[[obs_type]], sigma1s[obs_type])
                antigenic_map_short[,obs_type] <- create_cross_reactivity_vector(antigenic_map_melted[[obs_type]], sigma2s[obs_type])
            }

            y_new <- titre_data_fast(
                theta, theta_indices_unique, unique_obs_types,
                infection_history_mat, strain_isolation_times, infection_strain_indices,
                sample_times, type_data_start,obs_types,
                sample_data_start, titre_data_start,
                nrows_per_sample, measured_strain_indices, 
                antigenic_map_long,
                antigenic_map_short,
                antigenic_distances,
                mus, boosting_vec_indices,
                titre_before_infection
            )
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_par_tab]
                titre_shifts <- measurement_bias[expected_indices]
                y_new <- y_new + titre_shifts
            }
            y_new[overall_indices]
        }
    }
    f
}
