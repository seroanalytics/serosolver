#' Parallel tempering MCMC
#' 
#' Runs the MCMC chain for the serosolver model, but uses parallel tempering to run a number of chains at different temperatures, 
#'  swapping their states to improve mixing. NOTE that this can be quite slow!
#' @param par_tab as in \code{\link{run_MCMC}}, but a list of such data structures, one for each temperature, rather than a single data frame
#' @inheritParams run_MCMC
#' @details
#' The `mcmc_pars` argument has the following options in addition to those in \code{\link{run_MCMC}}:
#'  * temperature (a vector of temperatures to run the MCMC chains at)
#'  * parallel_tempering_iter how many iterations to run before swapping MCMC chains
#' @export
#' 
run_MCMC_pt <- function(par_tab,
                             titre_dat,
                             vaccination_histories=NULL,
                             vaccination_histories_mat=NULL,
                             antigenic_map=NULL,
                             strain_isolation_times=NULL,
                             mcmc_pars = c(),
                             mvr_pars = NULL,
                             start_inf_hist = NULL,
                             filename = "test",
                             CREATE_POSTERIOR_FUNC = create_posterior_func,
                             CREATE_PRIOR_FUNC = NULL,
                             version = 2,
                             mu_indices = NULL,
                             measurement_indices = NULL,
                             measurement_random_effects = FALSE,
                             proposal_ratios = NULL,
                             OPT_TUNING = 0.1,
                             solve_likelihood = TRUE,
                             n_alive = NULL,
                             continue_run = FALSE,
                             ...) {
    ## Error checks --------------------------------------
  if(!inherits(par_tab, "list")) {
    message(cat("List of par_tabs not given in parallel tempering sampler. Duplicating given table instead."))
    par_tab <- rep(list(par_tab), 10)
  }
  if(!inherits(start_inf_hist, "list")) {
    message(cat("List of start_inf_hist not given in parallel tempering sampler. Duplicating given table instead."))
    start_inf_hist <- rep(list(start_inf_hist), 10)
  }
  if(!("vac_flag" %in% par_tab[[1]]$names)) {
      for (i in 1:length(par_tab)) {
        par_tab[[i]] <- rbind(par_tab[[i]], list("vac_flag", 0, 1, 0.1, 0, 1,  0, 1, 0))
      }
  }
  lapply(par_tab, function(x) check_par_tab(x, TRUE, 2))
  

    ## Sort out MCMC parameters --------------------------------------
    ###################################################################
  mcmc_pars_used <- list(
    "iterations" = 20000, "popt" = 0.44, "popt_hist" = 0.44, "opt_freq" = 2000, "thin" = 1,
    "adaptive_period" = 5000,
    "save_block" = 100, "thin_hist" = 10, "hist_sample_prob" = 1, "switch_sample" = 2, "burnin" = 0,
    "inf_propn" = 0.5, "move_size" = 5, "hist_opt" = 0, "swap_propn" = 0.5,
    "hist_switch_prob" = 0, "year_swap_propn" = 1, "propose_from_prior"=TRUE,
    "temperature" = seq(1, 10, by = 1), "parallel_tempering_iter" = 5
  )

  mcmc_pars_used[names(mcmc_pars)] <- mcmc_pars
  cat(mcmc_pars_used[["temperature"]])

  mcmc_pars
  ## Extract MCMC parameters
  iterations <- mcmc_pars_used[["iterations"]] # How many iterations to run after adaptive period
  popt <- mcmc_pars_used[["popt"]] # Desired optimal acceptance rate
  popt_hist <- mcmc_pars_used[["popt_hist"]]
  opt_freq <- mcmc_pars_used[["opt_freq"]] # How often to adjust step size
  thin <- mcmc_pars_used[["thin"]] # Save only every nth iterations for theta sampling
  adaptive_period <- mcmc_pars_used[["adaptive_period"]] # How many iterations for adaptive period
  save_block <- mcmc_pars_used[["save_block"]] # How many post-thinning iterations to store before saving to disk
  hist_tab_thin <- mcmc_pars_used[["thin_hist"]] # Save only every nth iterations for infection history sampling
  hist_sample_prob <- mcmc_pars_used[["hist_sample_prob"]] # What proportion of infection histories to sample each step
  switch_sample <- mcmc_pars_used[["switch_sample"]] # Resample infection histories every n iterations
  burnin <- mcmc_pars_used[["burnin"]] # Run this many iterations before attempting adaptation. Idea is to reduce getting stuck in local maxima
  move_size <- mcmc_pars_used[["move_size"]] # Number of infections to move/remove/add in each proposal step
  inf_propn <- mcmc_pars_used[["inf_propn"]] # Number of infections to move/remove/add in each proposal step
  n_infs <- floor(length(antigenic_map$inf_times) * inf_propn)
  hist_opt <- mcmc_pars_used["hist_opt"] # Should infection history proposal step be adaptive?
  swap_propn <- mcmc_pars_used[["swap_propn"]] # If using gibbs, what proportion of proposals should be swap steps?
  hist_switch_prob <- mcmc_pars_used[["hist_switch_prob"]] # If using gibbs, what proportion of iterations should be swapping contents of two time periods?
  year_swap_propn <- mcmc_pars_used[["year_swap_propn"]] # If gibbs and swapping contents, what proportion of these time periods should be swapped?
  temperatures <- mcmc_pars_used[["temperature"]]
  parallel_tempering_iter <- mcmc_pars_used[["parallel_tempering_iter"]]
  propose_from_prior <- mcmc_pars_used[["propose_from_prior"]]
  ###################################################################
  # check if temps are monotonically increasing
  if (any(diff(temperatures) <= 0)) {
    stop("Temperature ladder is not monotonically increasing. Please change from: ", temperatures)
  }
  ###############
  ## Currently only doing this for gibbs version (2)
  ###############
  prior_on_total <- FALSE
  if (version == 1) { ## Lambda version
    stop("Prior version 1: There is not a parallel tempering sampler for this prior version yet")
  } else if (version == 2) { ## Gibbs version
    prop_print <- "Prior version 2: Using integrated FOI prior on infection history, with gibbs sampling of infections"
    hist_proposal <- 2
  } else if (version == 3) { ## Beta binomial version
    stop("Prior version 3: There is not a parallel tempering sampler for this prior version yet")
    hist_switch_prob <- 0
    hist_proposal <- 3
  } else if (version == 4) {
    stop("Prior version 4: There is not a parallel tempering sampler for this prior version yet")
    hist_proposal <- 2
    prior_on_total <- TRUE
  } else { ## By default, use phi version
    stop("Invalid version specified - 2 (beta on times) in the parallel tempering sampler")
  }
  if (!is.null(mvr_pars)) {
    stop("mvr_pars is not null - must be null as only univariate gibbs sampling available for parallel tempering sampler")
  }
  message(cat(prop_print, "\n"))

  mcmc_chain_file <- paste0(filename, "_chain.csv")
  infection_history_file <- paste0(filename, "_infection_histories.csv")
  mcmc_info_file <- paste0(filename, "_run_info.RData")
  check_data(titre_dat)

  ################################################################################################
  ###################### INITIALISE SECTION, ALWAYS NEEDS DOING ################################
  ################################################################################################

  par_tab_cold <- par_tab[[1]]
  par_names <- as.character(par_tab_cold$names)
  param_length <- nrow(par_tab_cold)
  unfixed_pars <- which(par_tab_cold$fixed == 0) 
  unfixed_par_length <- nrow(par_tab_cold[par_tab_cold$fixed == 0, ])
  lower_bounds <- par_tab_cold$lower_bound # Parameters cannot step below this
  upper_bounds <- par_tab_cold$upper_bound # Parameters cannot step above this
  alpha <- par_tab_cold[par_tab_cold$names == "alpha", "values"]
  beta <- par_tab_cold[par_tab_cold$names == "beta", "values"]
  
  ## If using phi terms, pull their indices out of the parameter table
  phi_indices <- NULL
  if ("phi" %in% par_names) {
    phi_indices <- which(par_tab$names == "phi")
  }

  if (!is.null(antigenic_map)) {
    strain_isolation_times <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
  } else {
    antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,"inf_times"=strain_isolation_times)
  } # How many strains are we testing against and what time did they circulate
  n_indiv <- length(unique(titre_dat$individual)) # How many individuals in the titre_dat?

  reset_par <- integer(param_length)
  reset_indiv <- integer(n_indiv)

  if (!is.null(titre_dat$DOB)) {
    DOBs <- unique(titre_dat[, c("individual", "DOB")])[, 2]
  } else {
    DOBs <- rep(min(strain_isolation_times), n_indiv)
  }
  age_mask <- create_age_mask(DOBs, strain_isolation_times)
  ## Create strain mask
  strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
  masks <- data.frame(cbind(age_mask, strain_mask))

  if (is.null(n_alive)) {
    n_alive <- get_n_alive_group(titre_dat, strain_isolation_times)
  }
  n_alive_tot <- sum(n_alive)

  if (switch_sample >= 1) {
    switch_sample_flag <- c(rep(TRUE, switch_sample), FALSE)
  } else {
    switch_sample_flag <- c(TRUE, rep(FALSE, floor(1 / switch_sample)))
  }

  switch_sample_flag_length <- length(switch_sample_flag)
  group_ids_vec <- unique(titre_dat[, c("individual", "group")])[, "group"] - 1
  n_groups <- length(unique(group_ids_vec))

  if (!is.null(mu_indices)) {
    prior_mu <- create_prior_mu(par_tab_cold)
  }
  if (measurement_random_effects) {
    prior_shifts <- create_prob_shifts(par_tab_cold)
  }

  extra_probabilities <- function(prior_pars, prior_infection_history) {
    names(prior_pars) <- par_names
    beta <- prior_pars["beta"]
    alpha <- prior_pars["alpha"]
    prior_probab <- 0
    ## If prior version 2 or 4
    if (hist_proposal == 2) {
      ## Prior version 4
      if (prior_on_total) {
        n_infections <- sum_infections_by_group(prior_infection_history, group_ids_vec, n_groups)
        n_infections_group <- rowSums(n_infections)
        prior_probab <- prior_probab + inf_mat_prior_total_group_cpp(
          n_infections_group,
          n_alive_tot, alpha, beta
        )
      } else {
          n_infections <- sum_infections_by_group(prior_infection_history, group_ids_vec, n_groups)
          if (any(n_infections > n_alive)) print("error")
          prior_probab <- prior_probab + inf_mat_prior_group_cpp(n_infections, n_alive, alpha, beta)
      }
    }
    if (!is.null(CREATE_PRIOR_FUNC)) prior_probab <- prior_probab + prior_func(prior_pars)
    if (!is.null(mu_indices)) prior_probab <- prior_probab + prior_mu(prior_pars)
    if (measurement_random_effects) prior_probab <- prior_probab + prior_shifts(prior_pars)
    prior_probab
  }

  # Define the two functions
  posterior_simp <- protect(CREATE_POSTERIOR_FUNC(
    par_tab_cold,
    titre_dat,
    vaccination_histories,
    vaccination_histories_mat,
    antigenic_map,
    strain_isolation_times,
    version = version,
    solve_likelihood,
    age_mask,
    measurement_indices_by_time = measurement_indices,
    mu_indices = mu_indices,
    n_alive = n_alive,
    function_type = 1,
    ...
  ))
  if (!is.null(CREATE_PRIOR_FUNC)) {
    prior_func <- CREATE_PRIOR_FUNC(par_tab_cold)
  }
  ## If using gibbs proposal on infection_history, create here
  proposal_gibbs <- protect(CREATE_POSTERIOR_FUNC(
    par_tab_cold,
    titre_dat,
    vaccination_histories,
    vaccination_histories_mat,
    antigenic_map,
    strain_isolation_times,
    version = version,
    solve_likelihood,
    age_mask,
    measurement_indices_by_time = measurement_indices,
    mu_indices = mu_indices,
    n_alive = n_alive,
    function_type = 2,
    ...
  ))

  if(!continue_run) {
    ################################################################################################
    ###################### INITIALISE SECTION, WITHOUT PRIOR SAMPLING ################################
    ################################################################################################
    #### ONLY NEED TO DO THIS IF NO RESMAPLING #######
    start_pars <- lapply(par_tab, function(x) x$values)
    steps_list_start <- lapply(par_tab, function(x) x$steps)
    par_tab_cold <- par_tab[[1]]
    current_pars <- par_tab_cold$values # Starting parameters, DONE
    steps <- par_tab_cold$steps # How far to step on unit scale to begin with?

    infection_history_swap_n <- infection_history_swap_accept <- 0
    switch_sample_i <- 1 # DONE

    if (is.null(mvr_pars)) {
      tempaccepted <- tempiter <- integer(param_length)
      reset <- integer(param_length)
      reset[] <- 0
    } else {
      tempaccepted <- tempiter <- reset <- 0
      cov_mat <- mvr_pars[[1]][unfixed_pars, unfixed_pars]
      steps <- mvr_pars[[2]]
      w <- mvr_pars[[3]]
    }


    histiter <- integer(n_indiv)
    histaccepted <- integer(n_indiv)
    histiter_add <- integer(n_indiv)
    histaccepted_add <- integer(n_indiv)
    histiter_move <- integer(n_indiv)
    histaccepted_move <- integer(n_indiv)
    
    overall_swap_proposals <- matrix(0, nrow = n_indiv, ncol = length(strain_isolation_times))
    overall_add_proposals <- matrix(0, nrow = n_indiv, ncol = length(strain_isolation_times))
    
    if(is.null(proposal_ratios)) {
        proposal_ratios <- rep(1, length(strain_isolation_times))
    }
    n_infs_vec <- rep(n_infs, n_indiv) # How many infection history moves to make with each proposal
    move_sizes <- rep(move_size, n_indiv) # How many years to move in smart proposal step

  # group_ids_vec <- unique(titre_dat[, c("individual", "group")])[, "group"] - 1
  # n_groups <- length(unique(group_ids_vec))

    ## Number of people that were born before each year and have had a sample taken since that year happened
  #  if (is.null(n_alive)) n_alive <- sapply(seq(1, length(strain_isolation_times)), function(x) nrow(masks[masks$age_mask <= x & masks$strain_mask >= x, ]))

    ## Create posterior calculating function

    ## If using random effects on mu, need to include hyperprior term on mu
    ## We can't do this in the main posterior function, because this term
    ## applies to the overall posterior whereas the main posterior function
    ## returns each individual's posterior

    ######################

    ## Setup initial conditions
    ## I think this needs to be a list
    infection_histories <- start_inf_hist[[1]]
    if (is.null(start_inf_hist[[1]])) {
        start_inf_hist <- lapply(1:length(temperatures), function(x) setup_infection_histories_titre(titre_dat, strain_isolation_times, space = 5, titre_cutoff = 3))
        infection_histories <- start_inf_hist[[1]]
    }

    check_inf_hist(titre_dat, strain_isolation_times, infection_histories)
    ## Initial likelihood
    tmp_posterior <- posterior_simp(current_pars, infection_histories)
    indiv_likelihoods <- tmp_posterior[[1]] / 1
    indiv_priors <- tmp_posterior[[2]]
      ## Initial total likelihoodss
    indiv_posteriors <- indiv_likelihoods + indiv_priors

      ## If needed for some proposal types per individual
    proposal_ratio <- rep(0, n_indiv)

    total_likelihood <- sum(indiv_likelihoods)
    total_prior_prob <- sum(indiv_priors) + extra_probabilities(current_pars, infection_histories)
    total_posterior <- total_likelihood + total_prior_prob
    
    message(cat("Starting posterior probability (coldest chain): ", total_posterior, "\n", sep = "\t"))
    message(cat("Starting likelihood (coldest chain): ", total_likelihood, "\n", sep = "\t"))
    message(cat("Starting prior prob (coldest chain): ", total_prior_prob, "\n", sep = "\t"))

    ## Initial indexing parameters
    sampno <- 2
    chain_index <- 1
    no_recorded <- 1
    par_i <- 1 # DONE
    i <- 1 # DONE
    i_prev <- 0
    log_prob <- 0
    cov_mat0 <- diag(unfixed_pars) # DONE


    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
    ## MCMC parameters (save these after sampling)
    ## Initialise values

    mcmc_list <- list(
      "i" = i, "par_i" = par_i, "current_pars" = current_pars,
      "infection_histories" = infection_histories,
      "indiv_likelihoods" = indiv_likelihoods, "total_likelihood" = total_likelihood,
      "total_prior_prob" = total_prior_prob, "total_posterior" = total_posterior,
      "tempaccepted" = tempaccepted, "tempiter" = tempiter,
      "steps" = steps, "temp" = 1,
      "cov_mat0" = cov_mat0,
      "histiter" = histiter,
      "histaccepted" = histaccepted,
      "histiter_add" = histiter_add,
      "histaccepted_add" = histaccepted_add,
      "histiter_move" = histiter_move,
      "histaccepted_move" = histaccepted_move,
      "overall_swap_proposals" = overall_swap_proposals,
      "overall_add_proposals" = overall_add_proposals,
      "proposal_ratios" = proposal_ratios,
      "propose_from_prior" = propose_from_prior,
      "infection_history_swap_accept" = infection_history_swap_accept,
      "infection_history_swap_n" = infection_history_swap_n,
      "switch_sample_i" = switch_sample_i,
      "move_sizes" = move_sizes,
      "n_infs_vec" = n_infs_vec
    )
    temperatures <- temperatures
    ## Replicate list for parallel tempering
    mcmc_list <- rep(list(mcmc_list), length(temperatures))
    ## Start values for parallel tempering

    mcmc_list <- Map(function(x, y) modifyList(x, list(current_pars = y)), mcmc_list, start_pars)
    mcmc_list <- Map(function(x, y) modifyList(x, list(steps = y)), mcmc_list, steps_list_start)
    mcmc_list <- Map(function(x, y) modifyList(x, list(temp = y)), mcmc_list, temperatures)
    mcmc_list <- Map(function(x, y) modifyList(x, list(infection_histories = y)), mcmc_list, start_inf_hist)
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

    # PRE ALLOCATE MEMORY
    ## Create empty chain to store every iteration for the adaptive period and burnin
    opt_chain <- matrix(nrow = burnin + adaptive_period, ncol = unfixed_par_length)
    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow = save_block, ncol = param_length + 4)

    ## Store posterior (called lnlike), likelihood ad prior
    chain_colnames <- c("sampno", par_names, "lnlike", "likelihood", "prior_prob")
    tmp_table <- array(dim = c(1, length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1, ] <- c(1, current_pars, total_posterior, total_likelihood, total_prior_prob)
    colnames(tmp_table) <- chain_colnames

    data.table::fwrite(as.data.frame(tmp_table),
      file = mcmc_chain_file,
      row.names = FALSE, col.names = TRUE, sep = ",", append = FALSE
    )
    
    save_infection_history_to_disk(infection_histories, infection_history_file, 1,
      append = FALSE, col_names = TRUE
    )
  } else {
    # IMPORT PREVIOUS RUN
    load(file = mcmc_info_file) # loads mcmc_info
    mcmc_list <- mcmc_info$mcmc_list
    temperatures <- mcmc_info$temperatures

    # PRE ALLOCATE MEMORY
    sampno <- mcmc_list[[1]]$i + 1
    i_prev <- mcmc_list[[1]]$i
    no_recorded <- 1
    chain_index <- 1

    ## Create empty chain to store every iteration for the adaptive period and burnin
    opt_chain <- matrix(nrow = burnin + adaptive_period, ncol = unfixed_par_length)
    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow = save_block, ncol = param_length + 4)
  }

  # Parallel tempering stuff
  offset <- 0
  potential_swaps <- swaps <- 0
  
  ## Main body of running MCMC
  for (i in (i_prev + 1):(iterations + adaptive_period + burnin + i_prev)) 
  {
        if (i %% save_block == 0) message(cat("Current iteration: ", i, "\n", sep = "\t"))
        for (jh in 1:length(mcmc_list)) mcmc_list[[jh]][["i"]] <- i
        for (k in 1:length(mcmc_list)) {
          ############################################
          ############################################
          ############################################
          ############### MCMC_PT_AUX START ###############
          ###########################################
          ############################################
          ############################################  
          mcmc_list_k <- mcmc_list[[k]] 
          par_i <- mcmc_list_k[["par_i"]]
          current_pars <- mcmc_list_k[["current_pars"]]
          infection_histories <- mcmc_list_k[["infection_histories"]]
          indiv_likelihoods <- mcmc_list_k[["indiv_likelihoods"]]
          total_likelihood <- mcmc_list_k[["total_likelihood"]]
          total_prior_prob <- mcmc_list_k[["total_prior_prob" ]]
          total_posterior <- mcmc_list_k[["total_posterior"]]
          tempaccepted <- mcmc_list_k[["tempaccepted"]]
          tempiter <- mcmc_list_k[["tempiter"]]
          steps <- mcmc_list_k[["steps" ]]
          temp <- mcmc_list_k[["temp"]]
          cov_mat0 <- mcmc_list_k[["cov_mat0"]]
          histiter <- mcmc_list_k[["histiter"]]
          histaccepted <- mcmc_list_k[["histaccepted"]]
          histiter_add <- mcmc_list_k[["histiter_add"]]
          histaccepted_add <- mcmc_list_k[["histaccepted_add"]]
          histiter_move <- mcmc_list_k[["histiter_move"]]
          histaccepted_move <- mcmc_list_k[["histaccepted_move"]]
          overall_swap_proposals <- mcmc_list_k[["overall_swap_proposals"]]
          overall_add_proposals <- mcmc_list_k[["overall_add_proposals"]]
          proposal_ratios <- mcmc_list_k[["proposal_ratios"]]
          propose_from_prior <- mcmc_list_k[["propose_from_prior"]]
          infection_history_swap_accept <- mcmc_list_k[["infection_history_swap_accept"]]
          infection_history_swap_n <- mcmc_list_k[["infection_history_swap_n"]]
          switch_sample_i <- mcmc_list_k[["switch_sample_i"]]
          move_sizes <- mcmc_list_k[["move_sizes"]]
          n_infs_vec <- mcmc_list_k[["n_infs_vec"]]

          ## Whether to swap entire year contents or not - only applies to gibbs sampling
          inf_swap_prob <- runif(1)

          ######################
          ## PROPOSALS
          ######################
          theta_sample <- switch_sample_flag[switch_sample_i]
          switch_sample_i <- switch_sample_i + 1
          if (switch_sample_i > switch_sample_flag_length) switch_sample_i <- 1
          
          ## A) EITHER SAMPLE A PARAMETER VALUE
          if (theta_sample) {
              ## If all pars are fixed
              if (length(unfixed_pars) == 0) {
                  proposal <- current_pars ## Set proposal to be current parameters
                  tempiter <- tempiter + 1
              } else {
                  ## If using univariate proposals
                  if (is.null(mvr_pars)) {
                  ## For each parameter (Gibbs)
                      j <- unfixed_pars[par_i]
                      par_i <- par_i + 1
                      if (par_i > unfixed_par_length) par_i <- 1
                      proposal <- univ_proposal(current_pars, lower_bounds, upper_bounds, steps, j)
                      tempiter[j] <- tempiter[j] + 1
                  ## If using multivariate proposals
                  } else {
                      proposal <- mvr_proposal(current_pars, unfixed_pars, steps * cov_mat,
                          steps * cov_mat0, FALSE,
                          beta = 0.05
                      )
                      tempiter <- tempiter + 1
                  }
              }
              ## Calculate new likelihood for these parameters
              tmp_new_posteriors <- posterior_simp(proposal, infection_histories)

              new_indiv_likelihoods <- tmp_new_posteriors[[1]] # For each individual
              new_indiv_priors <- tmp_new_posteriors[[2]]
              new_indiv_posteriors <- new_indiv_likelihoods + new_indiv_priors
              new_total_likelihood <- sum(new_indiv_likelihoods) # Total
              new_total_prior_prob <- sum(new_indiv_priors) +
                  extra_probabilities(proposal, infection_histories)
              new_total_posterior <- new_total_likelihood + new_total_prior_prob # Posterior
              ## B) RESAMPLE A INFECTION HISOTRY
          } else {

              ## Otherwise, resample infection history
              indiv_sub_sample <- sample(1:n_indiv, ceiling(hist_sample_prob * n_indiv))
              indiv_sub_sample <- indiv_sub_sample[order(indiv_sub_sample)]

              ## Generate random number 0-1 to decide whether to move an infection time, or add/remove one
              rand_ns <- runif(length(indiv_sub_sample))

              ## Need to temporarily store current parameters as new pars, as
              ## might change with phi swap step
              proposal <- current_pars
              names(proposal) <- par_names
              alpha <- proposal["alpha"]
              beta <- proposal["beta"]
              new_likelihoods_calculated <- FALSE

              if (inf_swap_prob > hist_switch_prob) {

                  prop_gibbs <- proposal_gibbs(
                      proposal,
                      infection_histories,
                      indiv_likelihoods,
                      indiv_sub_sample,
                      alpha, beta,
                      n_infs_vec, swap_propn, move_size,
                      histiter_add,
                      histaccepted_add,
                      histiter_move,
                      histaccepted_move,
                      overall_swap_proposals,
                      overall_add_proposals,
                      proposal_ratios,
                      temp,
                      propose_from_prior
                  )

                  histiter <- prop_gibbs$proposal_iter
                  histaccepted <- prop_gibbs$accepted_iter
                  new_indiv_likelihoods <- prop_gibbs$old_probs
                  new_infection_histories <- prop_gibbs$new_infection_history
                  new_likelihoods_calculated <- TRUE

                  overall_swap_proposals <- prop_gibbs$overall_swap_proposals
                  overall_add_proposals <- prop_gibbs$overall_add_proposals

              } else {
                  tmp <- inf_hist_swap(
                      infection_histories,
                      vaccination_histories_mat,
                      age_mask, strain_mask,
                      year_swap_propn, move_size
                  )
                  new_infection_histories <- tmp[[1]]
                  if (!identical(new_infection_histories, infection_histories)) {
                      infection_history_swap_n <- infection_history_swap_n + 1
                  }
              }
              ## The proposals are either a swap step or an add/remove step. Need to track which type was used for which individual,
              ## as we adapt the `step size` for these two update steps independently
              if (hist_proposal != 2 & inf_swap_prob > hist_switch_prob) {
                  move <- indiv_sub_sample[which(rand_ns < swap_propn)]
                  add <- indiv_sub_sample[which(rand_ns > swap_propn)]
                  histiter_add[add] <- histiter_add[add] + 1
                  histiter_move[move] <- histiter_move[move] + 1
                  histiter[indiv_sub_sample] <- histiter[indiv_sub_sample] + 1
              }
              ## Calculate new likelihood with these infection histories
              ## If we didn't calculate the new likelihoods above, then need to do so here
              #if (!new_likelihoods_calculated) {
                #  cat("new_infection_histories: ", new_infection_histories[1:100], "\n")
              new_post <- posterior_simp(proposal, new_infection_histories)
              new_indiv_likelihoods <- new_post[[1]]
              new_indiv_priors <- new_post[[2]]
            # }
              new_indiv_posteriors <- new_indiv_likelihoods + new_indiv_priors
              new_total_likelihood <- sum(new_indiv_likelihoods)
              new_total_prior_prob <- sum(new_indiv_priors) +
                  extra_probabilities(proposal, new_infection_histories)
              new_total_posterior <- new_total_likelihood + new_total_prior_prob
          }
          #############################
          ## METROPOLIS HASTINGS STEP
          #############################
          ## Check that all proposed parameters are in allowable range
          ## Skip if any parameters are outside of the allowable range
          log_prob <- new_total_posterior - total_posterior
          if (theta_sample) {
              if (!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)) {
                  log_prob <- min(log_prob, 0)

                  if (log(runif(1)) < log_prob / temp) {
                      if (!any(proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                          proposal[unfixed_pars] > upper_bounds[unfixed_pars])) {

                          ## Accept with probability 1 if better, or proportional to
                          ## difference if not
                          current_pars <- proposal
                          ## Store acceptances
                          ## If all parameters are fixed, then we 'accept'
                          if (length(unfixed_pars) == 0) {
                              tempaccepted <- tempaccepted + 1
                          } else {
                              if (is.null(mvr_pars)) {
                                  tempaccepted[j] <- tempaccepted[j] + 1
                              } else {
                                  tempaccepted <- tempaccepted + 1
                              }
                          }
                          indiv_likelihoods <- new_indiv_likelihoods
                          indiv_priors <- new_indiv_priors
                          indiv_posteriors <- new_indiv_posteriors

                          total_likelihood <- new_total_likelihood
                          total_prior_prob <- new_total_prior_prob
                          total_posterior <- new_total_posterior
                      }
                  }
              }
          } else {
              ## If normal proposal
              if (inf_swap_prob > hist_switch_prob) {
                  if (!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)) {
                      infection_histories <- new_infection_histories
                      indiv_likelihoods <- new_indiv_likelihoods
                      indiv_priors <- new_indiv_priors
                      indiv_posteriors <- new_indiv_posteriors
                      current_pars <- proposal
                      
                      total_likelihood <- new_total_likelihood
                      total_posterior <- new_total_posterior
                      total_prior_prob <- new_total_prior_prob
                  }
              } else {
                  if (!identical(new_infection_histories, infection_histories)) {
                      log_prob <- new_total_posterior - total_posterior
                      if (!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)) {
                          log_prob <- min(log_prob, 0)

                          if (log(runif(1)) < log_prob / temp) {
                              if (!any(proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                                  proposal[unfixed_pars] > upper_bounds[unfixed_pars])) {
                                  infection_history_swap_accept <- infection_history_swap_accept + 1
                                  infection_histories <- new_infection_histories
                                  current_pars <- proposal
                                  indiv_likelihoods <- new_indiv_likelihoods
                                  indiv_priors <- new_indiv_priors
                                  indiv_posteriors <- new_indiv_posteriors

                                  total_likelihood <- new_total_likelihood
                                  total_prior_prob <- new_total_prior_prob
                                  total_posterior <- new_total_posterior
                              }
                          }
                      }
                  }
              }
          }

          mcmc_list[[k]][["par_i"]] <- par_i
          mcmc_list[[k]][["current_pars"]] <- current_pars
          mcmc_list[[k]][["infection_histories"]] <- infection_histories
          mcmc_list[[k]][["indiv_likelihoods"]] <- indiv_likelihoods
          mcmc_list[[k]][["total_likelihood"]] <- total_likelihood
          mcmc_list[[k]][["total_prior_prob" ]] <- total_prior_prob
          mcmc_list[[k]][["total_posterior"]] <- total_posterior
          mcmc_list[[k]][["tempaccepted"]] <- tempaccepted
          mcmc_list[[k]][["tempiter"]] <- tempiter
          mcmc_list[[k]][["steps"]] <- steps
          mcmc_list[[k]][["temp"]] <- temp
          mcmc_list[[k]][["cov_mat0"]] <- cov_mat0
          mcmc_list[[k]][["histiter"]] <- histiter
          mcmc_list[[k]][["histaccepted"]] <- histaccepted 
          mcmc_list[[k]][["histiter_add"]] <- histiter_add
          mcmc_list[[k]][["histaccepted_add"]] <- histaccepted_add
          mcmc_list[[k]][["histiter_move"]] <- histiter_move
          mcmc_list[[k]][["histaccepted_move"]] <- histaccepted_move
          mcmc_list[[k]][["overall_swap_proposals"]] <- overall_swap_proposals
          mcmc_list[[k]][["overall_add_proposals"]] <- overall_add_proposals
          mcmc_list[[k]][["proposal_ratios"]] <- proposal_ratios
          mcmc_list[[k]][["propose_from_prior"]] <- propose_from_prior
          mcmc_list[[k]][["infection_history_swap_accept"]] <- infection_history_swap_accept
          mcmc_list[[k]][["infection_history_swap_n"]] <- infection_history_swap_n
          mcmc_list[[k]][["switch_sample_i"]] <- switch_sample_i
          mcmc_list[[k]][["move_sizes"]] <- move_sizes
          mcmc_list[[k]][["n_infs_vec"]] <- n_infs_vec

          ############################################
          ############################################
          ############################################
          ############### MCMC_PT_AUX END ###############
          ###########################################
          ############################################
          ############################################
        }
          ##############


        ## perform parallel tempering
        if (i %% parallel_tempering_iter == 0) {
            parallel_tempering_list <- parallel_tempering(mcmc_list, temperatures, offset)
            mcmc_list <- parallel_tempering_list$mcmc_list
            swaps <- swaps + parallel_tempering_list$swaps
            potential_swaps <- potential_swaps + 1
            offset <- 1 - offset
          }
          ##############################
          ## SAVE STEP
          ##############################
          ## If current iteration matches with recording frequency, store in the chain.
          ## If we are at the limit of the save block,
          ## save this block of chain to file and reset chain
          ## Save theta
        if (i %% thin == 0) {
            # use coldest chain
          current_pars <- mcmc_list[[1]][["current_pars"]]
          total_likelihood <- mcmc_list[[1]][["total_likelihood"]]
          total_prior_prob <- mcmc_list[[1]][["total_prior_prob"]]
          posterior <- mcmc_list[[1]][["total_posterior"]]
          save_chain[no_recorded, 1] <- sampno
          save_chain[no_recorded, 2:(ncol(save_chain) - 3)] <- current_pars
          save_chain[no_recorded, ncol(save_chain) - 2] <- posterior
          save_chain[no_recorded, ncol(save_chain) - 1] <- total_likelihood
          save_chain[no_recorded, ncol(save_chain)] <- total_prior_prob
          no_recorded <- no_recorded + 1
        }
        ## Save infection histories
        if (i %% hist_tab_thin == 0) {
          infection_histories <- mcmc_list[[1]][["infection_histories"]]
          save_infection_history_to_disk(infection_histories, infection_history_file, sampno)
        }

        ## USEFUL THINGS TO UPDATE
        ## Update steps
        scale_univariate <- function(steps, popt, pcur, unfixed_pars) {
          steps[unfixed_pars] <- vapply(
              unfixed_pars, function(x) scaletuning(steps[x], popt, pcur[x]),
              double(1)
          )
          steps
        }

        # pcur update
        tempaccepted_out <- lapply(mcmc_list, function(x) x[["tempaccepted"]])
        tempinter_out <- lapply(mcmc_list, function(x) x[["tempiter"]])

        pcur <- lapply(mcmc_list, function(x) x[["tempaccepted"]] / x[["tempiter"]])

        reset_infection_swap_pars <- function(mcmc_list, reset) {
          mcmc_list[["infection_history_swap_accept"]] <- mcmc_list[["infection_history_swap_n"]] <- reset
          mcmc_list
        }
        reset_acceptance <- function(mcmc_list, reset) {
            mcmc_list[["tempaccepted"]] <- mcmc_list[["tempiter"]] <- reset
            mcmc_list
        }
        reset_hist_pars <- function(mcmc_list, reset) {
            mcmc_list[["histaccepted"]] <- mcmc_list[["histiter"]] <- reset
            mcmc_list[["histaccepted_add"]] <- mcmc_list[["histiter_add"]] <- reset
            mcmc_list[["histaccepted_move"]] <- mcmc_list[["histiter_move"]] <- reset
            mcmc_list
        }

        ##############################
        ## ADAPTIVE PERIOD
        ##############################
        ## If within adaptive period, need to do some adapting!
        if ((i + i_prev) > (adaptive_period + burnin + i_prev) & (i %% opt_freq) == 0) {
          # Update step sizes
          steps_list <- lapply(mcmc_list, function(x) x[["steps"]])
          steps_list_updated <- Map(function(x, y) scale_univariate(x, popt, y, unfixed_pars), steps_list, pcur)
          mcmc_list <- Map(function(x, y) modifyList(x, list(steps = y)), mcmc_list, steps_list_updated)

          # Get coldest chain info
          infection_history_swap_accept <- mcmc_list[[1]]$infection_history_swap_accept
          infection_history_swap_n <- mcmc_list[[1]]$infection_history_swap_n
          histaccepted <- mcmc_list[[1]]$histaccepted
          histiter <- mcmc_list[[1]]$histiter
          histaccepted_add <- mcmc_list[[1]]$histaccepted_add
          histiter_add <- mcmc_list[[1]]$histiter_add
          histaccepted_move <- mcmc_list[[1]]$histaccepted_move
          histiter_move <- mcmc_list[[1]]$histiter_move

          pcur_hist <- histaccepted / histiter ## Overall
          pcur_hist_add <- histaccepted_add / histiter_add ## For adding
          pcur_hist_move <- histaccepted_move / histiter_move ## For adding
          pcur_hist_swap <- infection_history_swap_accept / infection_history_swap_n
          infection_histories <- mcmc_list[[1]][["infection_histories"]]

          # Print info during sampling
          message(cat("No infections: ", sum(infection_histories), "\n", sep = "\t"))
          message(cat("Pcur: ", signif(pcur[[1]][unfixed_pars], 3), "\n", sep = "\t"))
          message(cat("Step sizes: ", signif(steps_list_updated[[1]][unfixed_pars], 3), "\n", sep = "\t"))
          message(cat("Group inf hist swap pcur: ",
              signif(pcur_hist_swap, 3),"\n", 
              sep = "\t"
          ))

          swap_ratio <- swaps / potential_swaps
          message(cat("Swap ratio: ", swap_ratio, sep = "\t"))
          temperatures <- calibrate_temperatures(temperatures, swap_ratio)
          message(cat("Temperatures: ", temperatures, sep = "\t"))
          mcmc_list <- Map(function(x, y) modifyList(x, list(temp = y)), mcmc_list, temperatures)
          swaps <- potential_swaps <- 0

          message(cat("Pcur hist add: ", head(signif(pcur_hist_add, 3)),"\n",  sep = "\t"))
          message(cat("Pcur hist move: ", head(signif(pcur_hist_move, 3)), "\n", sep = "\t"))

          # Reset variables    
          mcmc_list <- lapply(mcmc_list, function(x) reset_acceptance(x, reset_par))
          mcmc_list <- lapply(mcmc_list, function(x) reset_infection_swap_pars(x, 0))
          mcmc_list <- lapply(mcmc_list, function(x) reset_hist_pars(x, reset_indiv))
        }
        if ((i + i_prev) > (burnin + i_prev) & (i  + i_prev) <= (adaptive_period + burnin + i_prev)) {
          opt_chain[chain_index, ] <- mcmc_list[[1]][["current_pars"]][unfixed_pars]
          if (chain_index %% opt_freq == 0) 
          {
            # Update step sizes
            steps_list <- lapply(mcmc_list, function(x) x[["steps"]])
            steps_list_updated <- Map(function(x, y) scale_univariate(x, popt, y, unfixed_pars), steps_list, pcur)
            mcmc_list <- Map(function(x, y) modifyList(x, list(steps = y)), mcmc_list, steps_list_updated)

            # Get coldest chain info
            histaccepted <- mcmc_list[[1]]$histaccepted
            histiter <- mcmc_list[[1]]$histiter
            histaccepted_add <- mcmc_list[[1]]$histaccepted_add
            histiter_add <- mcmc_list[[1]]$histiter_add
            histaccepted_move <- mcmc_list[[1]]$histaccepted_move
            histiter_move <- mcmc_list[[1]]$histiter_move
            infection_history_swap_accept <- mcmc_list[[1]]$infection_history_swap_accept
            infection_history_swap_n <- mcmc_list[[1]]$infection_history_swap_n

            pcur_hist <- histaccepted / histiter
            pcur_hist_add <- histaccepted_add / histiter_add
            pcur_hist_move <- histaccepted_move / histiter_move
            pcur_hist_swap <- infection_history_swap_accept / infection_history_swap_n

            ## Adaptive infection history proposal
            if (hist_opt == 1) {
              ## * Increase or decrease the number of infection history locations
              ## being changed to modify acceptance rate. If not accepting enough,
              ## reduce number. If accepting too many, increase number
              #for (k in 1:length(mcmc_list)) {
                k <- 1
                mcmc_list_k <- mcmc_list[[k]]
                move_sizes <- mcmc_list_k$move_sizes
                n_infs_vec <- mcmc_list_k$n_infs_vec
                n_infs_vec[which(pcur_hist_add < popt_hist * (1 - OPT_TUNING))] <-
                    n_infs_vec[which(pcur_hist_add < popt_hist * (1 - OPT_TUNING))] - 1
                n_infs_vec[which(pcur_hist_add >= popt_hist * (1 + OPT_TUNING))] <-
                    n_infs_vec[which(pcur_hist_add >= popt_hist * (1 + OPT_TUNING))] + 1
                n_infs_vec[n_infs_vec < 1] <- 1

                move_sizes[which(pcur_hist_move < popt_hist * (1 - OPT_TUNING))] <-
                    move_sizes[which(pcur_hist_move < popt_hist * (1 - OPT_TUNING))] - 1
                move_sizes[which(pcur_hist_move >= popt_hist * (1 + OPT_TUNING))] <-
                    move_sizes[which(pcur_hist_move >= popt_hist * (1 + OPT_TUNING))] + 1
                move_sizes[move_sizes < 1] <- 1

                for (ii in seq_along(n_infs_vec)) {
                    move_sizes[ii] <- min(move_sizes[ii], length(age_mask[ii]:strain_mask[ii]))
                    move_sizes[ii] <- min(move_sizes[ii], 10)
                    n_infs_vec[ii] <- min(n_infs_vec[ii], length(age_mask[ii]:strain_mask[ii]))
                }
                mcmc_list[[k]]$move_sizes <- move_sizes
                mcmc_list[[k]]$n_infs_vec <- n_infs_vec
            #  }
            }

            # Print info during sampling
            infection_histories <- mcmc_list[[1]][["infection_histories"]]
            ## Look at infection history proposal sizes
            message(cat("No infections: ", sum(infection_histories), "\n", sep = "\t"))
            message(cat("Pcur hist add: ", head(signif(pcur_hist_add, 3)), "\n", sep = "\t"))
            message(cat("No. infections sampled: ", head(n_infs_vec), "\n", sep = "\t"))
            message(cat("Pcur hist move: ", head(signif(pcur_hist_move, 3)), "\n", sep = "\t"))
            message(cat("Move sizes: ", head(move_sizes), "\n", sep = "\t"))
            message(cat("Pcur theta: ", signif(pcur[[1]][unfixed_pars], 3), "\n", sep = "\t"))
            message(cat("Step sizes: ", signif(steps_list_updated[[1]][unfixed_pars], 3), "\n", sep = "\t"))
            message(cat("Pcur group inf hist swap: ", signif(pcur_hist_swap, 3), "\n", sep = "\t"))
            message(cat("Group inf hist swap propn: ", year_swap_propn, "\n", sep = "\t"))
            ## If not accepting, send a warning
            if (all(pcur[[1]][unfixed_pars][!is.nan(pcur[[1]][unfixed_pars])] == 0)) {
                message("Warning: acceptance rates are 0. Might be an error with the theta proposal?\n")
            }
            # Reset variables    
            mcmc_list <- lapply(mcmc_list, function(x) reset_infection_swap_pars(x, 0))
            mcmc_list <- lapply(mcmc_list, function(x) reset_hist_pars(x, reset_indiv))
            mcmc_list <- lapply(mcmc_list, function(x) reset_acceptance(x, reset_par))
          }
          chain_index <- chain_index + 1
        }
        #######################
        ## HOUSEKEEPING
        #######################

        if (no_recorded == save_block) {
        data.table::fwrite(as.data.frame(save_chain[1:(no_recorded - 1), ]),
            file = mcmc_chain_file,
            col.names = FALSE, row.names = FALSE, sep = ",", append = TRUE
        )
        save_chain <- empty_save_chain
        no_recorded <- 1
        }
        sampno <- sampno + 1
  }
  ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
  ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
  ## rather than an array
  if (no_recorded > 2) {
      data.table::fwrite(as.data.frame(save_chain[1:(no_recorded - 1), ]),
                      file = mcmc_chain_file, row.names = FALSE,
                      col.names = FALSE, sep = ",", append = TRUE)
  }
  p_accept <- mcmc_list[[1]][["tempaccepted"]] / mcmc_list[[1]][["tempiter"]]
  p_accept <- p_accept[unfixed_pars]

  current_pars <- lapply(mcmc_list, function(x) x$current_pars)

  if (is.null(mvr_pars)) {
      cov_mat <- NULL
  }
  mcmc_info <- list(
    mcmc_list = mcmc_list,
    temperatures = temperatures
  )
  save(mcmc_info, file = mcmc_info_file)
  return(list(
      "chain_file" = mcmc_chain_file,
      "history_file" = infection_history_file,
      "cov_mat" = cov_mat,
      "current_pars" = current_pars,
      "steps" = steps,
      "temperatures" = temperatures,
      "overall_swap_proposals" = overall_swap_proposals,
      "overall_add_proposals" = overall_add_proposals
    )
  )
}


#' performs parallel tempering - Ada Yan
#'
#' @param mcmc_list a list of lists: values, log likelihood etc. of parallel MCMC chains
#' @param temperatures numeric vector: temperatures of chains
#' @param offset integer: 0 or 1. 0 = swap chains 1 <-> 2, 3 <-> 4...
#' 1 = swap chains 2<->3, 4<->5...
#'
#' @return a list of lists: values, log likelihood etc. of paralle chains after parallel tempering
#' @export
parallel_tempering <- function(mcmc_list, temperatures, offset) {
  recorded_swaps <- double(length(mcmc_list) - 1)

  ## extract current probabilities and log likelihoods

  all_likelihood <- vapply(mcmc_list, function(x) x$total_likelihood, double(1))
  all_prior_prob <- vapply(mcmc_list, function(x) x$total_prior_prob, double(1))
  all_posterior <- vapply(mcmc_list, function(x) x$total_posterior, double(1))
  all_likelihoods <- lapply(mcmc_list, function(x) x$indiv_likelihoods)
  all_steps <- lapply(mcmc_list, function(x) x$steps)

  all_current_pars <- lapply(mcmc_list, function(x) x$current_pars)
  all_infection_histories <- lapply(mcmc_list, function(x) x$infection_histories)

  ## decide which chains to swap
  if ((offset + 1) <= (length(mcmc_list) - 1)) {
    swap_ind <- seq(offset + 1, length(mcmc_list) - 1, by = 2)

    decide_if_swap <- function(x, y) {
      delta <- (1 / temperatures[y] - 1 / temperatures[x]) *
        (all_posterior[x] - all_posterior[y])
      runif(1) <= exp(delta)
    }

    swaps <- vapply(swap_ind, function(x) decide_if_swap(x, x + 1), logical(1))
    swap_ind <- swap_ind[swaps]

    ## perform swap
    perform_swap <- function(vec, swap_ind) {
      vec_new <- vec
      vec_new[swap_ind] <- vec[swap_ind + 1]
      vec_new[swap_ind + 1] <- vec[swap_ind]
      vec_new
    }
    all_current_pars <- perform_swap(all_current_pars, swap_ind)
    all_steps <- perform_swap(all_steps, swap_ind)
    all_likelihood <- perform_swap(all_likelihood, swap_ind)
    all_likelihoods <- perform_swap(all_likelihoods, swap_ind)
    all_prior_prob <- perform_swap(all_prior_prob, swap_ind)
    all_posterior <- perform_swap(all_posterior, swap_ind)
    all_infection_histories <- perform_swap(all_infection_histories, swap_ind)
    
    new_list <- Map(
      function(x, y, z, b, c, q, s) list(
          current_pars = x, infection_histories = y,
          indiv_likelihoods = z, total_likelihood = b,
          total_prior_prob = c, total_posterior = q,
          steps = s
        ),
      all_current_pars, all_infection_histories,
      all_likelihoods, all_likelihood,
      all_prior_prob, all_posterior, all_steps
    )

    mcmc_list <- Map(modifyList, mcmc_list, new_list)
    recorded_swaps[swap_ind] <- 1
  }
  list("swaps" = recorded_swaps, "mcmc_list" = mcmc_list)
}

#' calibrate temperatures for parallel chains - Ada Yan
#'
#' @param temperatures vector of length n: current temperatures of chains
#' @param swap_ratio vector of length n - 1: (proportion of accepted swaps
#' out of proposed swaps)/2
#' the factor of 2 arises because of the way the swaps are recorded, and how
#' we alternate between swapping 1<->2, 3<->4.... and 2<->3, 4<->5...
#' @return vector of length n: new temperatures of chains
#'
calibrate_temperatures <- function(temperatures, swap_ratio) {
  diff_temp <- diff(temperatures)
  ## find chains between which the swap ratio is too large
  too_large <- swap_ratio > .2 # note factor of 2 from main text -- see above
  ## find chains between which the swap ratio is too small
  too_small <- swap_ratio < .05
  ## adjust differences between temperatures accordingly
  diff_temp <- diff_temp * (too_large * 1.5 + too_small * .75 + (1 - too_large - too_small))

  diff_temp[diff_temp < 0.1] <- 0.1
  ## reconstruct temperatures from their differences
  
  cat(cumsum(c(temperatures[1], diff_temp)))
  cumsum(c(temperatures[1], diff_temp))
}