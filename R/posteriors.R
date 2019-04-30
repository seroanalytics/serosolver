#' Posterior function pointer
#'
#' Essential for the MCMC algorithm. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to pass the data set back and forth.
#' @param par_tab the parameter table controlling information such as bounds, initial values etc
#' @param titre_dat the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time)
#' @param antigenic_map a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_years
#' @param version which version of the posterior function to solve (corresponds mainly to the infection history prior). Mostly just left to 1, but there is one special case where this should be set to 4 for the gibbs sampler. This is only really used by \code{\link{run_MCMC}} to place the infection history prior on the total number of infections across all years and individuals when version = 4
#' @param solve_likelihood usually set to TRUE. If FALSE< does not solve the likelihood and instead just samples/solves based on the model prior
#' @param age_mask see \code{\link{create_age_mask}} - a vector with one entry for each individual specifying the first epoch of circulation in which an individual could have been exposed
#' @param measurement_indices_by_time if not NULL, then use these indices to specify which measurement bias parameter index corresponds to which time
#' @param mu_indices if not NULL, then use these indices to specify which boosting parameter index corresponds to which time
#' @param n_alive if not NULL, uses this as the number alive in a given year rather than calculating from the ages. This is needed if the number of alive individuals is known, but individual birth dates are not
#' @param function_type integer specifying which version of this function to use. Specify 1 to give a posterior solving function; 2 to give the gibbs sampler for infection history proposals; otherwise just solves the titre model and returns predicted titres. NOTE that this is not the same as the attack rate prior argument!
#' @param ... other arguments to pass to the posterior solving function
#' @return a single function pointer that takes only pars and infection_histories as unnamed arguments. This function goes on to return a vector of posterior values for each individual
#' @export
create_posterior_func <- function(par_tab,
                                  titre_dat,
                                  antigenic_map,
                                  version = 1,
                                  solve_likelihood = TRUE,
                                  age_mask = NULL,
                                  measurement_indices_by_time = NULL,
                                  mu_indices = NULL,
                                  n_alive = NULL,
                                  function_type = 1,
                                  ...) {
    ## Sort data in same way
    titre_dat <- titre_dat[order(titre_dat$individual, titre_dat$run, titre_dat$samples, titre_dat$virus), ]

    ## Isolate data table as vectors for speed
    titres <- titre_dat$titre

    ## Setup data vectors and extract
    setup_dat <- setup_titredat_for_posterior_func(titre_dat, antigenic_map, age_mask, n_alive)

    individuals <- setup_dat$individuals
    antigenic_map_melted <- setup_dat$antigenic_map_melted
    strain_isolation_times <- setup_dat$strain_isolation_times
    infection_strain_indices <- setup_dat$infection_strain_indices
    sample_times <- setup_dat$sample_times
    rows_per_indiv_in_samples <- setup_dat$rows_per_indiv_in_samples
    nrows_per_individual_in_data <- setup_dat$nrows_per_individual_in_data
    cum_nrows_per_individual_in_data <- setup_dat$cum_nrows_per_individual_in_data
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
    weights_indices <- which(par_tab$type == 4) ## For functional form version
  knot_indices <- which(par_tab$type == 5)
  mu_indices_par_tab <- which(par_tab$type == 6)
  #########################################################

  par_names_theta <- par_tab[theta_indices, "names"]

  ## Find which options are being used in advance for speed
  explicit_phi <- (length(phi_indices) > 0)
  spline_phi <- (length(knot_indices) > 0)
  use_measurement_bias <- (length(measurement_indices_par_tab) > 0) & !is.null(measurement_indices_by_time)
  titre_shifts <- NULL
  expected_indices <- NULL
  measurement_bias <- NULL
  use_strain_dependent <- (length(mu_indices) > 0) & !is.null(mu_indices)
  additional_arguments <- NULL

  if (use_measurement_bias) {
      expected_indices <- measurement_indices_by_time[match(titre_dat$virus, strain_isolation_times)]
  }

  if (use_strain_dependent) {
    additional_arguments <- list(
      "boosting_vec_indices" = mu_indices - 1,
      "mus" = rep(2, length(strain_isolation_times))
    )
  }

  ## Posterior calculating version
  if (function_type == 1) {
    ## If solving likelihood, calculate full likelihood/priors
    if (solve_likelihood) {
      ########################
      ## Function to return
      ########################
      f <- function(pars, infection_history_mat) {
        phis <- pars[phi_indices]
        theta <- pars[theta_indices]
        weights <- pars[weights_indices]
        knots <- pars[knot_indices]
        mus <- pars[mu_indices_par_tab]

        ## Add measurement shifts to titres
        if (use_measurement_bias) {
          measurement_bias <- pars[measurement_indices_par_tab]
          to_add <- measurement_bias[expected_indices]
        }

        ## Pass strain-dependent boosting down
        if (use_strain_dependent) {
          additional_arguments[["mus"]] <- mus
        }
        names(theta) <- par_names_theta
        
        ## Work out short and long term boosting cross reactivity - C++ function
        antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
        antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])

        ## Now pass to the C++ function
        y <- titre_data_group(
          theta, infection_history_mat, strain_isolation_times,
          infection_strain_indices, sample_times,
          rows_per_indiv_in_samples,
          cum_nrows_per_individual_in_data,
          nrows_per_blood_sample,
          measured_strain_indices,
          antigenic_map_long, antigenic_map_short,
          DOBs, additional_arguments
        )
        
        ## Calculate likelihoods and group per individual
        liks <- r_likelihood(y, titres, theta, expected_indices, measurement_bias)
        liks <- sum_buckets(liks, nrows_per_individual_in_data)

        ## If including explicit prior on FOI, add to individuals here
        if (explicit_phi) {
            liks <- liks + calc_phi_probs_indiv(phis, infection_history_mat, age_mask, strain_mask)
        }
        
        ## If using spline term for FOI, add here
        if (spline_phi) {
          liks <- liks + calc_phi_probs_monthly(
            phis, knots, weights,
            infection_history_mat, age_mask
          )
        }
        return(liks)
      }
      ## If not solving likelihood, just need priors
    } else {
      ## Function to return
      f <- function(pars, infection_history_mat) {
        phis <- pars[phi_indices]
        theta <- pars[theta_indices]
        weights <- pars[weights_indices]
        knots <- pars[knot_indices]
        mus <- pars[mu_indices_par_tab]

        liks <- rep(-100000, n_indiv)

        if (explicit_phi) {
          liks <- liks + calc_phi_probs_indiv(phis, infection_history_mat, age_mask, strain_mask)
        }
        if (spline_phi) {
          liks <- liks + calc_phi_probs_monthly(phis, knots, weights, infection_history_mat, age_mask)
        }

        return(liks)
      }
    }
    ## Gibbs proposal function
  } else if (function_type == 2) {
    ## Version 4 puts the prior on the total number of infections
    if (version == 4) {
      n_alive_total <- sum(n_alive)
    } else {
      n_alive_total <- -1
    }

    ## FUNCTION TO RETURN
    ## Gibbs proposal on infection histories
    f <- function(pars, infection_history_mat,
                      alpha, beta,
                      indiv_propn, n_years,
                      swap_propn = 0.5, swap_distance = 1, temp = 1) {
      
      phis <- pars[phi_indices]
      theta <- pars[theta_indices]
      weights <- pars[weights_indices]
      knots <- pars[knot_indices]
      mus <- pars[mu_indices_par_tab]
      names(theta) <- par_names_theta
      
      ## Measurement shifts
      if (use_measurement_bias) {
        measurement_bias <- pars[measurement_indices_par_tab]
        titre_shifts <- measurement_bias[expected_indices]
      }

      ## Strain dependent boosting, pass down
      if (use_strain_dependent) {
        additional_arguments[["mus"]] <- mus
      }
      names(theta) <- par_names_theta
      
      ## Work out short and long term boosting cross reactivity - C++ function
      antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
      antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])
      ## Now pass to the C++ function
      new_infection_history_mat <- infection_history_proposal_gibbs(
        theta,
        infection_history_mat,
        indiv_propn,
        n_years,
        age_mask,
        strain_mask,
        n_alive,
        swap_propn,
        swap_distance,
        alpha,
        beta,
        strain_isolation_times,
        infection_strain_indices,
        sample_times,
        rows_per_indiv_in_samples,
        cum_nrows_per_individual_in_data,
        nrows_per_blood_sample,
        measured_strain_indices,
        antigenic_map_long,
        antigenic_map_short,
        titres,
        titre_shifts,
        additional_arguments,
        DOBs,
        solve_likelihood,
        n_alive_total,
        temp
      )
      return(new_infection_history_mat)
    }
    ## Or model solving function
  } else {
    ## FUNCTION TO RETURN
    f <- function(pars, infection_history_mat) {
      phis <- pars[phi_indices]
      theta <- pars[theta_indices]
      weights <- pars[weights_indices]
      knots <- pars[knot_indices]
      mus <- pars[mu_indices_par_tab]

      if (use_strain_dependent) {
        additional_arguments[["mus"]] <- mus
      }

      names(theta) <- par_names_theta

      ## Work out short and long term boosting cross reactivity - C++ function
      antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
      antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])

      ## Now pass to the C++ function
      y <- titre_data_group(
        theta, infection_history_mat, strain_isolation_times,
        infection_strain_indices, sample_times,
        rows_per_indiv_in_samples,
        cum_nrows_per_individual_in_data,
        nrows_per_blood_sample,
        measured_strain_indices,
        antigenic_map_long, antigenic_map_short,
        DOBs, additional_arguments
      )

      ## Add measurement shifts
      if (use_measurement_bias) {
        measurement_bias <- pars[measurement_indices_par_tab]
        titre_shifts <- measurement_bias[expected_indices]
        y <- y + titre_shifts
      }

      return(y)
    }
  }
  f
}

#' Likelihood function given data
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

    ## Vectorise, calculate boundaries seperately
    liks <- numeric(length(expected))
    large_i <- data >= theta["MAX_TITRE"]
    small_i <- data < 1
    rest_i <- data >= 1 & data < theta["MAX_TITRE"]
                                        #large_i <- data > theta["MAX_TITRE"]
                                        #small_i <- data <= 0
                                        #rest_i <- data > 0 & data <= theta["MAX_TITRE"]

    liks[large_i] <- pnorm(theta["MAX_TITRE"], expected[large_i], theta["error"], lower.tail = FALSE, log.p = TRUE)
    liks[small_i] <- pnorm(1, expected[small_i], theta["error"], lower.tail = TRUE, log.p = TRUE)
    liks[rest_i] <- log(pnorm(data[rest_i] + 1, expected[rest_i], theta["error"], lower.tail = TRUE, log.p = FALSE) -
    pnorm(data[rest_i], expected[rest_i], theta["error"], lower.tail = TRUE, log.p = FALSE))
  return(liks)
}


#' Calculate FOI log probability
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @param phis a vector of FOIs
#' @param infection_history the matrix of infection histories
#' @param age_mask the age mask vector as returned by \code{\link{create_age_mask}}
#' @return a single log probability
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
#' @export
calc_phi_probs_indiv <- function(phis, infection_history, age_mask, strain_mask) {
  lik <- numeric(nrow(infection_history))
  for (i in 1:ncol(infection_history)) {
    lik <- lik + log(((phis[i]^infection_history[, i]) * (1 - phis[i])^(1 - infection_history[, i]))) * as.numeric(age_mask <= i) * as.numeric(strain_mask >= i)
  }
  lik
}

#' Calculate FOI from spline
#' 
#' Version of FOI prior that calculates a seasonal spline such that the FOI in a given time point (month) comes from this spline term
#' @inheritParams calc_phi_probs
#' @param foi vector of force of infections per year
#' @param knots vector of knots for the spline
#' @param theta vector of theta parameters for spline
#' @return a vector of log probabilities for each individual
#' @export
calc_phi_probs_monthly <- function(foi, knots, theta, infection_history, age_mask) {
  phis <- generate_phis(foi, knots, theta, length(foi), 12)
  lik <- numeric(nrow(infection_history))
  for (i in 1:ncol(infection_history)) {
    lik <- lik + infection_history[, i] * log(phis[i]) + (1 - infection_history[, i]) * log(phis[i]) + log(as.numeric(age_mask <= i)) + log(as.numeric(strain_mask >= i))
  }
  lik
}


#' FOR DEBUGGING
#'
#' Brute force implementation of calculating the explicit FOI
#' @export
calc_phi_probs_indiv_brute <- function(phis, infection_history, age_mask) {
  lik <- numeric(nrow(infection_history))
  for (j in 1:nrow(infection_history)) {
    lik[j] <- 0
    age <- age_mask[j]
    for (i in 1:ncol(infection_history)) {
      if (i >= age) {
        lik[j] <- lik[j] + log(phis[i]^infection_history[j, i] * (1 - phis[i])^(1 - infection_history[j, i]))
      }
    }
  }
  lik
}


#' Generate FOI phis
#'
#' Calculates the seasonal spline for \code{\link{calc_phi_probs_monthly}}
#' @inheritParams calc_phi_probs_monthly
#' @param n_years number of years to calculate
#' @param buckets number of buckets per year (12 for monthly, 1 for annual)
#' @param degree degree of the spline
#' @return a vector of FOIs for each time point
#' @export
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

#' Generates a spline for \code{\link{generate_phis}}
#' @export
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
#' @export
prob_mus <- function(mus, pars) {
  mu_mean <- pars["mu_mean"]
  mu_sd <- pars["mu_sd"]
  return(sum(dnorm(mus, mu_mean, mu_sd, log = TRUE)))
  # location <- log(mu_mean^2 / sqrt(mu_sd^2 + mu_mean^2))
  # shape <- sqrt(log(1 + (mu_sd^2/mu_mean^2)))
  # l_mean <- log(mu_mean) - (mu_sd^2)/2
  # p <- sum(dnorm(log(mus),mu_mean,mu_sd,log=TRUE))
  # p_mean <- 0.6
  # p_sd <- 0.5
  # p_mu <- log(p_mean/sqrt(1 + (p_sd/p_mean)^2))
  # p_sigma <- sqrt(log(1 + (p_sd/p_mean)^2))

  # p_lik <- log(p_sigma*2.506628) - 0.5*((mu_mean - p_mu)/p_sigma)^2

  # return(p+p_lik)
  # return(sum(log(dtruncnorm(mus, a=0,mean=mu_mean, sd=mu_sd))))
  # return(sum(dlnorm(mus, location, shape, log=TRUE)))
  # mean_log_y <- mean(log(mus))
  # sd_log_y <- sd(log(mus))
  # sigmaOfLogY <- dunif(mu_sd, 0.001*sd_log_y,1000*sd_log_y)
  # muOfLogY <- dnorm(mu_mean, mean_log_y, 1/(10*sd_log_y)^2)
  # return(sum(dlnorm(mus, mu_mean, 1/mu_sd^2, log=TRUE)) + sigmaOfLogY + muOfLogY)
  return(sum(dlnorm(mus, mu_mean, mu_sd, log = TRUE)))
}


#' Posterior function pointer - FAST VERSION
#'
#' Fast implementation of \code{\link{create_posterior_func}}. The old version is kept in for compatability, and for use where speed is not bottlenecked. Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to pass the data set back and forth.
#' @inheritParams create_posterior_func
#' @return a single function pointer that takes only pars and infection_histories as unnamed arguments. This function goes on to return a vector of posterior values for each individual
#' @export
create_posterior_func_fast <- function(par_tab,
                                       titre_dat,
                                       antigenic_map,
                                       version = 1,
                                       solve_likelihood = TRUE,
                                       age_mask = NULL,
                                       measurement_indices_by_time = NULL,
                                       mu_indices = NULL,
                                       n_alive = NULL,
                                       function_type = 1,
                                       ...) {


    ## Seperate out initial readings and repeat readings - we only
    ## want to solve the model once for each unique indiv/sample/virus year tested
    titre_dat_unique <- titre_dat[titre_dat$run == 1, ]
    titre_dat_repeats <- titre_dat[titre_dat$run != 1, ]
    tmp <- row.match(titre_dat_repeats[, c("individual", "samples", "virus")],
                     titre_dat_unique[, c("individual", "samples", "virus")])
    titre_dat_repeats$index <- tmp

    overall_indices <- row.match(titre_dat[,c("individual", "samples", "virus")],
                                 titre_dat_unique[, c("individual", "samples", "virus")])
    
    ## Setup data vectors and extract
    setup_dat <- setup_titredat_for_posterior_func(titre_dat_unique, antigenic_map, age_mask, n_alive)

    individuals <- setup_dat$individuals
    antigenic_map_melted <- setup_dat$antigenic_map_melted
    strain_isolation_times <- setup_dat$strain_isolation_times
    infection_strain_indices <- setup_dat$infection_strain_indices
    sample_times <- setup_dat$sample_times
    rows_per_indiv_in_samples <- setup_dat$rows_per_indiv_in_samples
    nrows_per_individual_in_data <- setup_dat$nrows_per_individual_in_data
    cum_nrows_per_individual_in_data <- setup_dat$cum_nrows_per_individual_in_data
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
    weights_indices <- which(par_tab$type == 4) ## For functional form version
    knot_indices <- which(par_tab$type == 5)
    mu_indices_par_tab <- which(par_tab$type == 6)
#########################################################
    
    ## Some additional setup for the repeat data
    nrows_per_individual_in_data_repeats <- NULL
    for (individual in unique(individuals)) {
        nrows_per_individual_in_data_repeats <- c(nrows_per_individual_in_data_repeats, nrow(titre_dat_repeats[titre_dat_repeats$individual == individual, ]))
    }
    cum_nrows_per_individual_in_data_repeats <- cumsum(c(0, nrows_per_individual_in_data_repeats))

    titres_unique <- titre_dat_unique$titre
    titres_repeats <- titre_dat_repeats$titre
    repeat_indices <- titre_dat_repeats$index
    repeat_indices_cpp <- repeat_indices - 1

    par_names_theta <- par_tab[theta_indices, "names"]

    ## Find which options are being used in advance for speed
    explicit_phi <- (length(phi_indices) > 0)
    spline_phi <- (length(knot_indices) > 0)
    use_measurement_bias <- (length(measurement_indices_par_tab) > 0) & !is.null(measurement_indices_by_time)
    titre_shifts <- c(0)
    expected_indices <- NULL
    measurement_bias <- NULL
    use_strain_dependent <- (length(mu_indices) > 0) & !is.null(mu_indices)
    additional_arguments <- NULL

    if (use_measurement_bias) {
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

  if (function_type == 1) {
      f <- function(pars, infection_history_mat) {
          theta <- pars[theta_indices]
          names(theta) <- par_names_theta

          if (use_strain_dependent) {
              mus <- pars[mu_indices_par_tab]
          }
          
          antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma1"])
          antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["sigma2"])

          y_new <- titre_data_fast(
              theta, infection_history_mat, strain_isolation_times, infection_strain_indices,
              sample_times, rows_per_indiv_in_samples, cum_nrows_per_individual_in_data,
              nrows_per_blood_sample, measured_strain_indices, antigenic_map_long,
              antigenic_map_short, mus, boosting_vec_indices
          )

          if (use_measurement_bias) {
              measurement_bias <- pars[measurement_indices_par_tab]
              titre_shifts <- measurement_bias[expected_indices]
              y_new <- y_new + titre_shifts
          }
          #return(list(theta, titres_unique, y_new, titres_repeats, repeat_indices, nrows_per_individual_in_data, nrows_per_individual_in_data_repeats))
          if(solve_likelihood){
              ## Calculate likelihood for unique titres and repeat data
              liks <- likelihood_func_fast(theta, titres_unique, y_new)
              liks_repeats <- likelihood_func_fast(theta, titres_repeats, y_new[repeat_indices])
              
              ## Sum these for each individual
              liks <- sum_buckets(liks, nrows_per_individual_in_data) +
                  sum_buckets(liks_repeats, nrows_per_individual_in_data_repeats)
              ## If including explicit prior on FOI, add to individuals here

              if (explicit_phi) {
                  phis <- pars[phi_indices]
                  liks <- liks + calc_phi_probs_indiv(phis, infection_history_mat, age_mask, strain_mask)
              }
              
              ## If using spline term for FOI, add here
              if (spline_phi) {
                  weights <- pars[weights_indices]
                  knots <- pars[knot_indices]
                  liks <- liks + calc_phi_probs_monthly(
                                     phis, knots, weights,
                                     infection_history_mat, age_mask
                                 )
              }
          } else {
              liks <- rep(-100000, n_indiv)
          }
          liks
      }
  } else if (function_type == 2) {
      f <- function(pars, infection_history_mat,
                    probs, sampled_indivs,
                    alpha, beta,
                  n_infs, swap_propn, swap_dist,
                  temp) {
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
        
        ## Now pass to the C++ function
        res <- infection_history_proposal_gibbs_fast(
            theta,
            infection_history_mat,
            probs,
            sampled_indivs,
            n_infs,
            age_mask,
            strain_mask,
            n_alive,
            swap_propn,
            swap_dist,
            alpha,
            beta,
            strain_isolation_times,
            infection_strain_indices,
            sample_times,
            rows_per_indiv_in_samples,
            cum_nrows_per_individual_in_data,
            cum_nrows_per_individual_in_data_repeats,
            nrows_per_blood_sample,
            measured_strain_indices,
            antigenic_map_long,
            antigenic_map_short,
            titres_unique,
            titres_repeats,
            repeat_indices_cpp,
            titre_shifts,
            mus,
            boosting_vec_indices,
            temp,
            solve_likelihood
        )
        return(res)
    }
  } else {
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
              antigenic_map_short, mus, boosting_vec_indices
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

#' Calculate likelihood of infection history matrix and sensitivity given inf_data
#'
#' @param inf_dat a matrix of the same dimension as infection_histories, 0 (not infected) and 1 (infected)
#' @param infection_histories infection history matrix
#' @param pars the current parameter values
#' @param par_tab  parameter table controlling information such as bounds, initial values etc
#' @return log likelihood value
#' @export
inf_likelihood<-function(inf_dat,infection_histories,pars,par_tab){
  
  theta_indices <- which(par_tab$type %in% c(0, 1))
  par_names_theta <- par_tab[theta_indices, "names"]
  theta <- pars[theta_indices]
  names(theta) <- par_names_theta
  
  #extract sensitivity parameter from parTab 
  rho <- theta["delta"]
  
  #of the infected in X
  inf_vec<-infection_histories[which(infection_histories==1)]
  #what's the match with Y
  obs_vec<-inf_dat[which(inf_dat==1)]
  
  #calculate LOG likleihood
  lik<-sum(dbinom(x=obs_vec,size=inf_vec,prob=rho,log=T),na.rm=T)
  
  return(lik)
}
