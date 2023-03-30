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
  rest_i <- data >= theta["MIN_TITRE"] & data <= theta["MAX_TITRE"]
  
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
#' @param titre_dat the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time). See \code{\link{example_titre_dat}}
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
                                  ...) {
    #browser()
    check_par_tab(par_tab, TRUE, version)
    if (!("group" %in% colnames(titre_dat))) {
        titre_dat$group <- 1
    }
    check_data(titre_dat)
   
    if (!is.null(antigenic_map)) {
      strain_isolation_times <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
    } else {
      antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,"inf_times"=strain_isolation_times)
    }
    
    ## Seperate out initial readings and repeat readings - we only
    ## want to solve the model once for each unique indiv/sample/virus year tested
    titre_dat_unique <- titre_dat[titre_dat$run == 1, ]
    ## Observations from repeats
    titre_dat_repeats <- titre_dat[titre_dat$run != 1, ]
    ## Find which entry in titre_dat_unique each titre_dat_repeats entry should correspond to
    tmp <- row.match(
        titre_dat_repeats[, c("individual", "samples", "virus")],
        titre_dat_unique[, c("individual", "samples", "virus")]
    )
    titre_dat_repeats$index <- tmp


    ## Which entries in the overall titre_dat matrix does each entry in titre_dat_unique correspond to?
    overall_indices <- row.match(
        titre_dat[, c("individual", "samples", "virus")],
        titre_dat_unique[, c("individual", "samples", "virus")]
    )
    ## Setup data vectors and extract
    setup_dat <- setup_titredat_for_posterior_func(
        titre_dat_unique, antigenic_map, 
        strain_isolation_times,
        age_mask, n_alive
    )

    individuals <- setup_dat$individuals
    n_groups <- length(unique(titre_dat$group))
    antigenic_map_melted <- setup_dat$antigenic_map_melted
    antigenic_distances <- c(melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")]))
    strain_isolation_times <- setup_dat$strain_isolation_times
    infection_strain_indices <- setup_dat$infection_strain_indices
    sample_times <- setup_dat$sample_times
    rows_per_indiv_in_samples <- setup_dat$rows_per_indiv_in_samples
    nrows_per_individual_in_data <- setup_dat$nrows_per_individual_in_data
    cum_nrows_per_individual_in_data <- setup_dat$cum_nrows_per_individual_in_data
    group_id_vec <- setup_dat$group_id_vec

    nrows_per_blood_sample <- setup_dat$nrows_per_blood_sample
    measured_strain_indices <- setup_dat$measured_strain_indices
    n_alive <- setup_dat$n_alive
    age_mask <- setup_dat$age_mask
    strain_mask <- setup_dat$strain_mask
    n_indiv <- setup_dat$n_indiv
    DOBs <- setup_dat$DOBs


#########################################################
    ## Extract parameter type indices from par_tab, to split up
    ## similar parameters in model solving functions
    option_indices <- which(par_tab$type == 0)
    theta_indices <- which(par_tab$type %in% c(0, 1))
    phi_indices <- which(par_tab$type == 2)
    measurement_indices_par_tab <- which(par_tab$type == 3)
    weights_indices <- which(par_tab$type == 4) ## For functional form version of FOI
    knot_indices <- which(par_tab$type == 5) ## For functional form version of FOI
    mu_indices_par_tab <- which(par_tab$type == 6)
#########################################################

    ## Some additional setup for the repeat data
    nrows_per_individual_in_data_repeats <- NULL
    nrows_per_individual_in_data_repeats <- plyr::ddply(titre_dat, .(individual),
                                                        function(x) nrow(x[x$run != 1,]))$V1
    cum_nrows_per_individual_in_data_repeats <- cumsum(c(0, nrows_per_individual_in_data_repeats))

    titres_unique <- titre_dat_unique$titre
    titres_repeats <- titre_dat_repeats$titre
    repeat_indices <- titre_dat_repeats$index
    repeat_indices_cpp <- repeat_indices - 1

    par_names_theta <- par_tab[theta_indices, "names"]

    ## Find which options are being used in advance for speed
    explicit_phi <- (length(phi_indices) > 0) ## Explicit prob of infection term
    spline_phi <- (length(knot_indices) > 0)
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

    if (function_type == 1) {
        message(cat("Creating posterior solving function...\n"))
        f <- function(pars, infection_history_mat) {
            theta <- pars[theta_indices]
            names(theta) <- par_names_theta

            if (use_strain_dependent) {
                mus <- pars[mu_indices_par_tab]
            }

            antigenic_map_long <- create_cross_reactivity_vector(
                antigenic_map_melted,
                theta["sigma1"]
            )
            antigenic_map_short <- create_cross_reactivity_vector(
                antigenic_map_melted,
                theta["sigma2"]
            )

            ## Calculate titres for measured data
            y_new <- titre_data_fast(
                theta, infection_history_mat, strain_isolation_times, infection_strain_indices,
                sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data,
                nrows_per_blood_sample, measured_strain_indices,
                antigenic_map_long,
                antigenic_map_short,
                antigenic_distances,
                mus, boosting_vec_indices
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
                liks <- likelihood_func_fast(theta, titres_unique, y_new)
                liks <- sum_buckets(liks, nrows_per_individual_in_data)
                if (repeat_data_exist) {
                    liks_repeats <- likelihood_func_fast(theta, titres_repeats, y_new[repeat_indices])
                    liks <- liks + sum_buckets(liks_repeats, nrows_per_individual_in_data_repeats)
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
                nrows_per_blood_sample,
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
                solve_likelihood
            )
            return(res)
        }
    } else {
        message(cat("Creating model solving function...\n"))
        ## Final version is just the model solving function
        f <- function(pars, infection_history_mat) {
            theta <- pars[theta_indices]
            names(theta) <- par_names_theta

            ## Pass strain-dependent boosting down
            if (use_strain_dependent) {
                mus <- pars[mu_indices_par_tab]
            }

            antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
            antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])

            y_new <- titre_data_fast(
                theta, infection_history_mat, strain_isolation_times, infection_strain_indices,
                sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data,
                nrows_per_blood_sample, measured_strain_indices, antigenic_map_long,
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
