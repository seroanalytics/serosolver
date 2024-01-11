#' Adaptive Metropolis-within-Gibbs/Metropolis Hastings Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. The algorithm then uses \code{\link{univ_proposal}} or \code{\link{mvr_proposal}} to explore parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename. This version of the algorithm is also designed to explore posterior densities for infection histories. See the package vignettes for examples. 
#' @param par_tab The parameter table controlling information such as bounds, initial values etc. See \code{\link{example_par_tab}}
#' @param antibody_data The data frame of titre data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); titre (integer of titre value against the given virus at that sampling time); run (integer giving the repeated number of this titre); DOB (integer giving date of birth matching time units used in model). See \code{\link{example_antibody_data}}
#' @param antigenic_map (optional) A data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. See \code{\link{example_antigenic_map}}
#' @param possible_exposure_times (optional) If no antigenic map is specified, this argument gives the vector of times at which individuals can be infected
#' @param mcmc_pars Named vector named vector with parameters for the MCMC procedure. See details
#' @param mvr_pars Leave NULL to use univariate proposals. Otherwise, a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @param start_inf_hist Infection history matrix to start MCMC at. Can be left NULL. See \code{\link{example_inf_hist}}
#' @param filename The full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param CREATE_POSTERIOR_FUNC Pointer to posterior function used to calculate a likelihood. This will probably be \code{\link{create_posterior_func}}
#' @param CREATE_PRIOR_FUNC User function of prior for model parameters. Should take parameter values only
#' @param prior_version which infection history assumption prior_version to use? See \code{\link{describe_priors}} for options. Can be 1, 2, 3 or 4
#' @param measurement_indices optional NULL. For measurement bias function. Vector of indices of length equal to number of circulation times. For each year, gives the index of parameters named "rho" that correspond to each time period
#' @param measurement_random_effects optional FALSE. Boolean indicating if measurement bias is a random effects term. If TRUE adds a component to the posterior calculation that calculates the probability of a set of measurement shifts "rho", given a mean and standard deviation
#' @param proposal_ratios optional NULL. Can set the relative sampling weights of the infection state times. Should be an integer vector of length matching nrow(antigenic_map). Otherwise, leave as NULL for uniform sampling.
#' @param OPT_TUNING Constant describing the amount of leeway when adapting the proposals steps to reach a desired acceptance rate (ie. does not change step size if within OPT_TUNING of the specified acceptance rate)
#' @param temp Temperature term for parallel tempering, raises likelihood to this value. Just used for testing at this point
#' @param solve_likelihood if FALSE, returns only the prior and does not solve the likelihood. Use this if you wish to sample directly from the prior
#' @param n_alive if not NULL, uses this as the number alive for the infection history prior, rather than calculating the number alive based on antibody_data
#' @param ... Other arguments to pass to CREATE_POSTERIOR_FUNC
#' @return A list with: 1) relative file path at which the MCMC chain is saved as a .csv file; 2) relative file path at which the infection history chain is saved as a .csv file; 3) the last used covariance matrix if mvr_pars != NULL; 4) the last used scale/step size (if multivariate proposals) or vector of step sizes (if univariate proposals)
#' @details
#' The `mcmc_pars` argument has the following options:
#'  * iterations (number of post adaptive period iterations to run)
#'  * adaptive_period (for this many iterations, change proposal step size adaptively every `opt_freq` iterations)
#'  * opt_freq (adapt proposal step size every opt_freq iterations)
#'  * thin (save every n iterations)
#'  * thin_hist (infection history thinning)
#'  * hist_sample_prob (proportion of individuals resampled at each iteration)
#'  * save_block (after this many iterations (post thinning), save to disk)
#'  * popt (desired acceptance rate for parameters)
#'  * popt_hist (desired acceptance rate for infection histories)
#'  * switch_sample (ratio of infection history samples to theta samples ie. switch_sample = 2 means sample inf hist twice for every theta sample, switch_sample = 0.5 means sample theta twice for every inf hist sample)
#'  * inf_propn (proportion of infection times to resample for each individual at each iteration)
#'  * move_size (number of infection years/months to move when performing infection history swap step)
#'  * hist_opt (if 1, performs adaptive infection history proposals. If 0, retains the starting infection history proposal parameters. This should ONLY be turned on for prior version 2)
#'  * swap_propn (if using gibbs sampling of infection histories, what proportion of proposals should be swap steps)
#'  * hist_switch_prob (proportion of infection history proposal steps to swap year_swap_propn of two time periods' contents)
#'  * year_swap_propn (when swapping contents of two time points, what proportion of individuals should have their contents swapped)
#' @md
#' @seealso \url{https://github.com/jameshay218/lazymcmc}
#' @family mcmc
#' @examples
#' \dontrun{
#' data(example_antibody_data)
#' data(example_par_tab)
#' data(example_antigenic_map)
#' res <- run_MCMC(example_par_tab[example_par_tab$names != "phi",], example_antibody_data, example_antigenic_map, prior_version=2)
#' }
#' @export
run_MCMC <- function(par_tab,
                     antibody_data,
                     antigenic_map=NULL,
                     possible_exposure_times=NULL,
                     mcmc_pars = c(),
                     mvr_pars = NULL,
                     start_inf_hist = NULL,
                     filename = "test",
                     CREATE_POSTERIOR_FUNC = create_posterior_func,
                     CREATE_PRIOR_FUNC = NULL,
                     prior_version = 1,
                     measurement_indices = NULL,
                     measurement_random_effects = FALSE,
                     proposal_ratios = NULL,
                     OPT_TUNING = 0.1,
                     temp = 1,
                     solve_likelihood = TRUE,
                     n_alive = NULL,
                     ...) {
  ## Error checks --------------------------------------
  if (!is.null(antigenic_map)) {
    possible_exposure_times <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
  } else {
    antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,"inf_times"=possible_exposure_times)
  }
  
  par_tab <- check_par_tab(par_tab, TRUE,possible_exposure_times=possible_exposure_times, version=prior_version)
  ## Sort out MCMC parameters --------------------------------------
  ###################################################################
  mcmc_pars_used <- c(
    "iterations" = 50000, "popt" = 0.44, "popt_hist" = 0.44, "opt_freq" = 2000, "thin" = 1,
    "adaptive_period" = 10000,
    "save_block" = 100, "thin_hist" = 10, "hist_sample_prob" = 0.5, "switch_sample" = 2, "burnin" = 0,
    "inf_propn" = 0.5, "move_size" = 3, "hist_opt" = 0, "swap_propn" = 0.5,
    "hist_switch_prob" = 0, "year_swap_propn" = 1, "propose_from_prior"=TRUE
  )
    mcmc_pars_used[names(mcmc_pars)] <- mcmc_pars

    ## Extract MCMC parameters
    iterations <- mcmc_pars_used["iterations"] # How many iterations to run after adaptive period
    popt <- mcmc_pars_used["popt"] # Desired optimal acceptance rate
    popt_hist <- mcmc_pars_used["popt_hist"]
    opt_freq <- mcmc_pars_used["opt_freq"] # How often to adjust step size
    thin <- mcmc_pars_used["thin"] # Save only every nth iterations for theta sampling
    adaptive_period <- mcmc_pars_used["adaptive_period"] # How many iterations for adaptive period
    save_block <- mcmc_pars_used["save_block"] # How many post-thinning iterations to store before saving to disk
    hist_tab_thin <- mcmc_pars_used["thin_hist"] # Save only every nth iterations for infection history sampling
    hist_sample_prob <- mcmc_pars_used["hist_sample_prob"] # What proportion of infection histories to sample each step
    switch_sample <- mcmc_pars_used["switch_sample"] # Resample infection histories every n iterations
    burnin <- mcmc_pars_used["burnin"] # Run this many iterations before attempting adaptation. Idea is to reduce getting stuck in local maxima
    move_size <- mcmc_pars_used["move_size"] # Number of infections to move/remove/add in each proposal step
    inf_propn <- mcmc_pars_used["inf_propn"] # Number of infections to move/remove/add in each proposal step
    n_infs <- floor(length(antigenic_map$inf_times) * inf_propn)
    hist_opt <- mcmc_pars_used["hist_opt"] # Should infection history proposal step be adaptive?
    swap_propn <- mcmc_pars_used["swap_propn"] # If using gibbs, what proportion of proposals should be swap steps?
    hist_switch_prob <- mcmc_pars_used["hist_switch_prob"] # If using gibbs, what proportion of iterations should be swapping contents of two time periods?
    year_swap_propn <- mcmc_pars_used["year_swap_propn"] # If gibbs and swapping contents, what proportion of these time periods should be swapped?
    propose_from_prior <- mcmc_pars_used["propose_from_prior"]
  ###################################################################

  ## Sort out which prior_version to run --------------------------------------
  prior_on_total <- FALSE
  if (prior_version == 1) { ## Lambda prior_version
    prop_print <- "Prior version 1: Using phi prior on infection history, with symmetric proposal probabilities"
    hist_proposal <- 1
  } else if (prior_version == 2) { ## Gibbs prior_version
    prop_print <- "Prior version 2: Using integrated FOI prior on infection history, with gibbs sampling of infections"
    hist_proposal <- 2
  } else if (prior_version == 3) { ## Beta binomial prior_version
    prop_print <- "Prior version 3: Using beta prior on total number of infections for an individual, with proposals from this prior"
    hist_switch_prob <- 0
    hist_proposal <- 3
  } else if (prior_version == 4) {
    prop_print <- "Prior version 4: Using beta prior on total number of infections across all years and all individuals, with gibbs sampling of infections"
    hist_proposal <- 2
    prior_on_total <- TRUE
  } else { ## By default, use phi prior_version
    stop("Invalid version specified - must be 1 (phi), 2 (beta on times), 3 (beta on individual) or 4 (beta on overall)")
  }
  message(cat(prop_print, "\n"))

  ## Extract parameter settings
  par_names <- as.character(par_tab$names) # Parameter names

  param_length <- nrow(par_tab)
  unfixed_pars <- which(par_tab$fixed == 0) # Indices of free parameters
  unfixed_par_length <- nrow(par_tab[par_tab$fixed == 0, ]) # How many free parameters?
  current_pars <- par_tab$values # Starting parameters

  ## Parameter constraints
  lower_bounds <- par_tab$lower_bound # Parameters cannot step below this
  upper_bounds <- par_tab$upper_bound # Parameters cannot step above this
  steps <- par_tab$steps # How far to step on unit scale to begin with? "steps" will be added above by check_par_tab

  ## If using phi terms, pull their indices out of the parameter table
  phi_indices <- NULL
  if ("phi" %in% par_names) {
    phi_indices <- which(par_tab$names == "phi")
  }

  infection_model_prior_shape1 <- par_tab[par_tab$names == "infection_model_prior_shape1", "values"]
  infection_model_prior_shape2 <- par_tab[par_tab$names == "infection_model_prior_shape2", "values"]

  ## To store acceptance rate of entire time period infection history swaps
  infection_history_swap_n <- infection_history_swap_accept <- 0
  ## Arrays to store acceptance rates
  ## If univariate proposals, store vector of acceptances
  if (is.null(mvr_pars)) {
    tempaccepted <- tempiter <- integer(param_length)
    reset <- integer(param_length)
    reset[] <- 0
  } else {
    ## If multivariate proposals, aggregate to one acceptance rate.
    ## Also extract covariance matrix, scale of proposal steps, and how
    ## much weighting to give to previous covariance matrix upon adaptive update
    tempaccepted <- tempiter <- reset <- 0
    cov_mat <- mvr_pars[[1]][unfixed_pars, unfixed_pars]
    steps <- mvr_pars[[2]]
    w <- mvr_pars[[3]]
  }
  ## Setup MCMC chain file with correct column names
  mcmc_chain_file <- paste0(filename, "_chain.csv")
  infection_history_file <- paste0(filename, "_infection_histories.csv")


  ###############
  ## Extract antibody_data parameters
  ##############
  ## Check the antibody_data input
  antibody_data <- check_data(antibody_data)
  
  n_indiv <- length(unique(antibody_data$individual)) # How many individuals in the antibody_data?

   
###################
    ## Housekeeping for infection history chain
###################
    histiter <- integer(n_indiv)
    histaccepted <- integer(n_indiv)
    histiter_add <- integer(n_indiv)
    histaccepted_add <- integer(n_indiv)
    histiter_move <- integer(n_indiv)
    histaccepted_move <- integer(n_indiv)

    overall_swap_proposals <- matrix(0,nrow=n_indiv,ncol=length(possible_exposure_times))
    overall_add_proposals <- matrix(0,nrow=n_indiv,ncol=length(possible_exposure_times))

    ## Scaling of infection history time proposal sample probs    
    if(is.null(proposal_ratios)){
        proposal_ratios <- rep(1, length(possible_exposure_times))
    }
        
    ##histadd_overall <- integer(n_indiv)
    ##histmove_overall <- integer(n_indiv)

    n_infs_vec <- rep(n_infs, n_indiv) # How many infection history moves to make with each proposal
    move_sizes <- rep(move_size, n_indiv) # How many years to move in smart proposal step
###############
    ## Create age mask
    ## -----------------------
  ## Note that DOBs for all groups must be from same reference point
  ## -----------------------
  ###############
  if (!is.null(antibody_data$birth)) {
    DOBs <- unique(antibody_data[, c("individual", "birth")])[, 2]
  } else {
    DOBs <- rep(min(possible_exposure_times), n_indiv)
  }
  age_mask <- create_age_mask(DOBs, possible_exposure_times)
  ## Create strain mask
  sample_mask <- create_sample_mask(antibody_data, possible_exposure_times)
  masks <- data.frame(cbind(age_mask, sample_mask))

  group_ids_vec <- unique(antibody_data[, c("individual", "population_group")])[, "population_group"] - 1
  n_groups <- length(unique(group_ids_vec))
  ## Number of people that were born before each year and have had a sample taken since that year happened

  if (is.null(n_alive)) {
    n_alive <- get_n_alive_group(antibody_data, possible_exposure_times)
  }
  ## Create posterior calculating function
  posterior_simp <- protect(CREATE_POSTERIOR_FUNC(par_tab,
    antibody_data,
    antigenic_map,
    possible_exposure_times,
    prior_version=prior_version,
    solve_likelihood,
    age_mask,
    measurement_indices_by_time = measurement_indices,
    n_alive = n_alive,
    function_type = 1,
    ...
  ))

  if (!is.null(CREATE_PRIOR_FUNC)) {
    prior_func <- CREATE_PRIOR_FUNC(par_tab)
  }

  ## If using gibbs proposal on infection_history, create here
  if (hist_proposal == 2) {
    proposal_gibbs <- protect(CREATE_POSTERIOR_FUNC(par_tab,
      antibody_data,
      antigenic_map,
      possible_exposure_times,
      prior_version=prior_version,
      solve_likelihood,
      age_mask,
      measurement_indices_by_time = measurement_indices,
      n_alive = n_alive,
      function_type = 2,
      ...
    ))
  }
    
    if (measurement_random_effects) {
        prior_shifts <- create_prob_shifts(par_tab)
    }
######################

    ## Setup initial conditions
    infection_histories <- start_inf_hist

    if (is.null(start_inf_hist)) {
        infection_histories <- setup_infection_histories_antibody_level(antibody_data, possible_exposure_times, space = 5, antibody_cutoff = 3)
    }

  
    check_inf_hist(antibody_data, possible_exposure_times, infection_histories)

    ## Initial likelihoods and individual priors
    tmp_posterior <- posterior_simp(current_pars, infection_histories)
    indiv_likelihoods <- tmp_posterior[[1]] / temp
    indiv_priors <- tmp_posterior[[2]]
    ## Initial total likelihoods
    indiv_posteriors <- indiv_likelihoods + indiv_priors

    ## If needed for some proposal types per individual
    proposal_ratio <- rep(0, n_indiv)
    n_alive_tot <- rowSums(n_alive)
    ## Create closure to add extra prior probabilities, to avoid re-typing later
    extra_probabilities <- function(prior_pars, prior_infection_history) {
        names(prior_pars) <- par_names
        infection_model_prior_shape1 <- prior_pars["infection_model_prior_shape1"]
        infection_model_prior_shape2 <- prior_pars["infection_model_prior_shape2"]
        prior_probab <- 0

        ## If prior prior_version 2 or 4
        if (hist_proposal == 2) {
          ## Prior prior_version 4
          if (prior_on_total) {
            n_infections <- sum_infections_by_group(prior_infection_history, group_ids_vec, n_groups)
            n_infections_group <- rowSums(n_infections)
            prior_probab <- prior_probab + inf_mat_prior_total_group_cpp(
              n_infections_group,
              n_alive_tot, infection_model_prior_shape1, infection_model_prior_shape2
            )
          } else {
              n_infections <- sum_infections_by_group(prior_infection_history, group_ids_vec, n_groups)
              if (any(n_infections > n_alive)){
                 #message("error -- more infections than there are individuals alive")
                 prior_probab <- -Inf
              } else {
                prior_probab <- prior_probab + 
                  inf_mat_prior_group_cpp(n_infections, n_alive, infection_model_prior_shape1, infection_model_prior_shape2)
              }
          }
        }
        if (!is.null(CREATE_PRIOR_FUNC)) prior_probab <- prior_probab + prior_func(prior_pars)
        if (measurement_random_effects) prior_probab <- prior_probab + prior_shifts(prior_pars)
        prior_probab
    }
    ## Initial total prior prob
    total_prior_prob <- sum(indiv_priors) + extra_probabilities(current_pars,infection_histories)
    total_likelihood <- sum(indiv_likelihoods)
    ## Initial posterior prob
    total_posterior <- total_likelihood + total_prior_prob

    if(!is.finite(total_prior_prob)) stop(paste("Error: starting prior probability is not finite."))
    if(!is.finite(total_likelihood)) stop(paste("Error: starting likelihood is not finite."))
    if(!is.finite(total_posterior)) stop(paste("Error: starting posterior probability is not finite."))
    
    message(cat("Starting posterior probability: ", total_posterior, "\n", sep = "\t"))
    message(cat("Starting likelihood : ", total_likelihood, "\n", sep = "\t"))
    message(cat("Starting prior prob: ", total_prior_prob, "\n", sep = "\t"))
###############

  ####################
  ## PRE ALLOCATE MEMORY
  ####################
  ## Create empty chain to store every iteration for the adaptive period and burnin
  opt_chain <- matrix(nrow = burnin + adaptive_period, ncol = unfixed_par_length)

  ## Create empty chain to store "save_block" iterations at a time
  save_chain <- empty_save_chain <- matrix(nrow = save_block, ncol = param_length + 4)

  ## Set up initial csv file
  ## Store posterior (called posterior_prob), likelihood ad prior
  chain_colnames <- c("sampno", par_names, "posterior_prob", "likelihood", "prior_prob")
  tmp_table <- array(dim = c(1, length(chain_colnames)))
  tmp_table <- as.data.frame(tmp_table)
  tmp_table[1, ] <- c(1, current_pars, total_posterior, total_likelihood, total_prior_prob)
  colnames(tmp_table) <- chain_colnames

  ## Write starting conditions to file
  data.table::fwrite(as.data.frame(tmp_table),
    file = mcmc_chain_file,
    row.names = FALSE, col.names = TRUE, sep = ",", append = FALSE
  )

  save_infection_history_to_disk(infection_histories, infection_history_file, 1,
    append = FALSE, col_names = TRUE
  )
  ## Initial indexing parameters
  no_recorded <- 1
  sampno <- 2
  par_i <- 1
  chain_index <- 1

  cov_mat0 <- diag(unfixed_pars)
  #####################
  ## MCMC ALGORITHM
  #####################
  log_prob <- 0
  
  ## If switch_sample is 0, only ever sample theta. If Inf, only sample inf_hist
  ## Track whether we should be sampling theta or inf hist
  if (switch_sample == 0){
      switch_sample_flag <- TRUE
  } else if (switch_sample >= 1) {
    switch_sample_flag <- c(rep(TRUE, switch_sample), FALSE)
  } else if(is.infinite(switch_sample)) {
      switch_sample_flag <- FALSE
  } else  {
    switch_sample_flag <- c(TRUE, rep(FALSE, floor(1 / switch_sample)))
  }
  
  switch_sample_i <- 1
  switch_sample_flag_length <- length(switch_sample_flag)

  for (i in 1:(iterations + adaptive_period + burnin)) {
    ## Whether to swap entire year contents or not - only applies to gibbs sampling
    inf_swap_prob <- runif(1)
    if (i %% save_block == 0) message(cat("Current iteration: ", i, "\n", sep = "\t"))
  
    ######################
    ## PROPOSALS
    ######################
    ## Keep track of switching between theta and inf hist
    theta_sample <- switch_sample_flag[switch_sample_i]
    switch_sample_i <- switch_sample_i + 1
    if (switch_sample_i > switch_sample_flag_length) switch_sample_i <- 1

    ## If updating theta
    if (theta_sample) {
      ## If all pars are fixed
      if (length(unfixed_pars) == 0) {
        proposal <- current_pars ## Set proposal to be current parameters
        tempiter <- tempiter + 1
      } else { ## Else propose parameters
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
      new_indiv_likelihoods <- tmp_new_posteriors[[1]] / temp # For each individual
      new_indiv_priors <- tmp_new_posteriors[[2]]
      new_indiv_posteriors <- new_indiv_likelihoods + new_indiv_priors
      new_total_likelihood <- sum(new_indiv_likelihoods) # Total
      new_total_prior_prob <- sum(new_indiv_priors) + extra_probabilities(proposal, infection_histories)
      new_total_posterior <- new_total_likelihood + new_total_prior_prob # Posterior

        ## Otherwise, resample infection history
    } else {
        ## Choose a random subset of individuals to update
        indiv_sub_sample <- sample(1:n_indiv, ceiling(hist_sample_prob * n_indiv))
        indiv_sub_sample <- indiv_sub_sample[order(indiv_sub_sample)]

        ## Generate random number 0-1 to decide whether to move an infection time, or add/remove one
        rand_ns <- runif(length(indiv_sub_sample))

        ## Need to temporarily store current parameters as new pars, as
        ## might change with phi swap step
        proposal <- current_pars
        names(proposal) <- par_names
        infection_model_prior_shape1 <- proposal["infection_model_prior_shape1"]
        infection_model_prior_shape2 <- proposal["infection_model_prior_shape2"]
        new_likelihoods_calculated <- FALSE ## Flag if we calculate the new likelihoods earlier than anticipated
        ## Which infection history proposal to use?
        ## Explicit phis on infection histories
        if (hist_proposal == 1) {
            ## Either swap entire contents or propose new infection history matrix
            if (inf_swap_prob > hist_switch_prob) {
                new_infection_histories <- infection_history_symmetric(
                    infection_histories,
                    indiv_sub_sample,
                    age_mask, sample_mask, move_sizes,
                    n_infs_vec, rand_ns, swap_propn
                )
            } else {
                tmp_hist_switch <- inf_hist_swap_phi(
                    infection_histories, proposal[phi_indices],
                    age_mask, sample_mask, year_swap_propn,
                    move_size, n_alive
                )
                new_infection_histories <- tmp_hist_switch[[1]]
                proposal[phi_indices] <- tmp_hist_switch[[2]]
                infection_history_swap_n <- infection_history_swap_n + 1
            }
            ## Gibbs sampler prior_version, integrate out phi
        } else if (hist_proposal == 2) {
            ## Swap entire contents or propose new
            if (inf_swap_prob > hist_switch_prob) {
                prop_gibbs <- proposal_gibbs(
                    proposal,
                    infection_histories,
                    indiv_likelihoods,
                    indiv_sub_sample,
                    infection_model_prior_shape1, infection_model_prior_shape2,
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
                new_indiv_likelihoods <- prop_gibbs$likelihoods_pre_proposal_tmp
                new_infection_histories <- prop_gibbs$new_infection_history     
                
                new_likelihoods_calculated <- TRUE

                overall_swap_proposals <- prop_gibbs$overall_swap_proposals
                overall_add_proposals <- prop_gibbs$overall_add_proposals
            } else {
                tmp <- inf_hist_swap(
                    infection_histories,
                    age_mask, sample_mask,
                    year_swap_propn, move_size
                )
                new_infection_histories <- tmp[[1]]
                if (!identical(new_infection_histories, infection_histories)) {
                    infection_history_swap_n <- infection_history_swap_n + 1
                }
                new_likelihoods_calculated <- FALSE
            }
            ## Beta binomial on per individual total infections
        } else if (hist_proposal == 3) {
            new_infection_histories <- inf_hist_prop_prior_v3(
                infection_histories, indiv_sub_sample,
                age_mask, sample_mask,
                move_sizes, n_infs_vec,
                infection_model_prior_shape1, infection_model_prior_shape2, rand_ns, swap_propn
            )
            ## Otherwise, default to symmetric proposals (hopefully this is never called)
        } else {
            new_infection_histories <- infection_history_symmetric(
                infection_histories, indiv_sub_sample,
                age_mask, sample_mask, move_sizes, n_infs_vec,
                rand_ns, swap_propn
            )
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
        if (!new_likelihoods_calculated) {
            new_post <- posterior_simp(proposal, new_infection_histories)
            new_indiv_likelihoods <- new_post[[1]] / temp
            new_indiv_priors <- new_post[[2]]
        }
        
        new_indiv_posteriors <- new_indiv_likelihoods + new_indiv_priors
        new_total_likelihood <- sum(new_indiv_likelihoods)
        new_total_prior_prob <- sum(new_indiv_priors) + extra_probabilities(proposal, new_infection_histories)
        new_total_posterior <- new_total_likelihood + new_total_prior_prob
    }
    #############################
    ## METROPOLIS HASTINGS STEP
    #############################
    ## Check that all proposed parameters are in allowable range
    ## Skip if any parameters are outside of the allowable range
    log_prob <- new_total_posterior - total_posterior
    #if(!is.finite(new_total_prior_prob)) browser()
    if (theta_sample) {
      if (!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)) {
        log_prob <- min(log_prob, 0)
        if (log(runif(1)) < log_prob) {
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
            ## MH step for each individual
            ## For prior prior_version 1 and 3, need to explicitly do acceptances
            if (hist_proposal != 2) {
                log_probs <- (new_indiv_posteriors[indiv_sub_sample] -
                              indiv_posteriors[indiv_sub_sample]) +
                    proposal_ratio[indiv_sub_sample]
                log_probs[log_probs > 0] <- 0
                x <- which(log(runif(length(indiv_sub_sample))) < log_probs)
                change_i <- indiv_sub_sample[x]

                infection_histories[change_i, ] <- new_infection_histories[change_i, ]

                indiv_likelihoods[change_i] <- new_indiv_likelihoods[change_i]
                indiv_priors[change_i] <- new_indiv_priors[change_i]
                indiv_posteriors[change_i] <- new_indiv_posteriors[change_i]

                total_likelihood <- sum(indiv_likelihoods)
                total_prior_prob <- sum(indiv_priors) +
                    extra_probabilities(current_pars, infection_histories)
                total_posterior <- total_likelihood + total_prior_prob

                ## Record acceptances for each add or move step
                add <- intersect(add, change_i)
                move <- intersect(move, change_i)

                histaccepted_add[add] <- histaccepted_add[add] + 1
                histaccepted_move[move] <- histaccepted_move[move] + 1
                histaccepted[change_i] <- histaccepted[change_i] + 1

                proposal_ratio <- rep(0, n_indiv)
                ## For prior prior_versions 2 and 4, have already done the acceptance rates
            } else {
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
            }
            ## Otherwise, doing the alternative swapping function
      } else {
        if (!identical(new_infection_histories, infection_histories)) {
          log_prob <- new_total_posterior - total_posterior
          if (!is.na(log_prob) & !is.nan(log_prob) & is.finite(log_prob)) {
            log_prob <- min(log_prob, 0)
            if (log(runif(1)) < log_prob) {
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


    ##############################
    ## SAVE STEP
    ##############################
    ## If current iteration matches with recording frequency, store in the chain.
    ## If we are at the limit of the save block,
    ## save this block of chain to file and reset chain

    ## Save theta
    if (i %% thin == 0) {
      save_chain[no_recorded, 1] <- sampno
      save_chain[no_recorded, 2:(ncol(save_chain) - 3)] <- current_pars
      save_chain[no_recorded, ncol(save_chain) - 2] <- total_posterior
      save_chain[no_recorded, ncol(save_chain) - 1] <- total_likelihood
      save_chain[no_recorded, ncol(save_chain)] <- total_prior_prob
      no_recorded <- no_recorded + 1
    }

    ## Save infection histories
    if (i %% hist_tab_thin == 0) {
      save_infection_history_to_disk(infection_histories, infection_history_file, sampno)
    }

    ##############################
    ## ADAPTIVE PERIOD
    ##############################
    ## If within adaptive period, need to do some adapting!
    if (i > (adaptive_period + burnin) & i %% opt_freq == 0) {
      pcur <- tempaccepted / tempiter ## get current acceptance rate
      message(cat("Pcur: ", signif(pcur, 3), "\n", sep = "\t"))
      message(cat("Step sizes: ", signif(steps, 3), "\n", sep = "\t"))
      message(cat("Group inf hist swap pcur: ",
        signif(infection_history_swap_accept / infection_history_swap_n, 3),"\n", 
        sep = "\t"
      ))
      infection_history_swap_accept <- infection_history_swap_n <- 0
      tempaccepted <- tempiter <- reset
      ## Have a look at the acceptance rates for infection histories
      ## if(hist_proposal != 2){
      pcur_hist <- histaccepted / histiter ## Overall
      pcur_hist_add <- histaccepted_add / histiter_add ## For adding
      pcur_hist_move <- histaccepted_move / histiter_move ## For adding

      message(cat("Pcur hist add: ", head(signif(pcur_hist_add, 3)),"\n",  sep = "\t"))
      message(cat("Pcur hist move: ", head(signif(pcur_hist_move, 3)), "\n", sep = "\t"))

      ##histadd_overall <- histadd_overall + histiter_add
      ##histmove_overall <- histmove_overall + histiter_move
      
      histiter <- integer(n_indiv)
      histaccepted <- integer(n_indiv)
      histiter_add <- integer(n_indiv)
      histaccepted_add <- integer(n_indiv)
      histiter_move <- integer(n_indiv)
      histaccepted_move <- integer(n_indiv)
      ## }
    }
    if (i > burnin & i <= (adaptive_period + burnin)) {
      ## Current acceptance rate
      pcur <- tempaccepted / tempiter
      ## Save each step
      opt_chain[chain_index, ] <- current_pars[unfixed_pars]
      ## If in an adaptive step
      if (chain_index %% opt_freq == 0) {
        ## If using univariate proposals
        if (is.null(mvr_pars)) {
          ## For each non fixed parameter, scale the step size
          for (x in unfixed_pars) steps[x] <- scaletuning(steps[x], popt, pcur[x])
        } else {
          if (chain_index > OPT_TUNING * adaptive_period & chain_index < adaptive_period) {
            old_cov_mat <- cov_mat
            ## Creates a new covariance matrix, but weights it with the old one
            cov_mat <- cov(opt_chain[1:chain_index, ])
            cov_mat <- w * cov_mat + (1 - w) * old_cov_mat
          }
          ## Scale tuning for last 80% of the adaptive period
          if (chain_index > (0.2) * adaptive_period) {
              steps <- scaletuning(steps, popt, pcur)
          }
        }
          pcur_hist <- histaccepted / histiter
          
          ##histadd_overall <- histadd_overall + histiter_add
          ##histmove_overall <- histmove_overall + histiter_move
          
          pcur_hist_add <- histaccepted_add / histiter_add
          pcur_hist_move <- histaccepted_move / histiter_move
          pcur_hist_swap <- infection_history_swap_accept / infection_history_swap_n
          
          infection_history_swap_accept <- infection_history_swap_n <- 0
          histiter <- integer(n_indiv)
          histaccepted <- integer(n_indiv)
          histiter_add <- integer(n_indiv)
          histaccepted_add <- integer(n_indiv)
          histiter_move <- integer(n_indiv)
          histaccepted_move <- integer(n_indiv)
          if (hist_opt == 1) {
              ## If adaptive infection history proposal
              ## Increase or decrease the number of infection history locations
              ## being changed to modify acceptance rate. If not accepting enough,
              ## reduce number. If accepting too many, increase number
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
                  move_sizes[ii] <- min(move_sizes[ii], length(age_mask[ii]:sample_mask[ii]))
                  move_sizes[ii] <- min(move_sizes[ii], 10)
                  n_infs_vec[ii] <- min(n_infs_vec[ii], length(age_mask[ii]:sample_mask[ii]))
              }
          }
          ## Look at infection history proposal sizes
          message(cat("Pcur hist add: ", head(signif(pcur_hist_add, 3)), "\n", sep = "\t"))
          message(cat("No. infections sampled: ", head(n_infs_vec), "\n", sep = "\t"))
          message(cat("Pcur hist move: ", head(signif(pcur_hist_move, 3)), "\n", sep = "\t"))
          message(cat("Move sizes: ", head(move_sizes), "\n", sep = "\t"))
          message(cat("Pcur theta: ", signif(pcur, 3), "\n", sep = "\t"))
          message(cat("Step sizes: ", signif(steps, 3), "\n", sep = "\t"))
          message(cat("Pcur group inf hist swap: ", signif(pcur_hist_swap, 3), "\n", sep = "\t"))
          message(cat("Group inf hist swap propn: ", year_swap_propn, "\n", sep = "\t"))
          pcur_hist <- histaccepted / histiter ## Overall
          ## If not accepting, send a warning
          if (all(pcur[!is.nan(pcur)] == 0)) {
              message("Warning: acceptance rates are 0. Might be an error with the theta proposal?\n")
          }

          tempaccepted <- tempiter <- reset
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
                           file = mcmc_chain_file, row.names = FALSE, col.names = FALSE,
                           sep = ",", append = TRUE
                           )
    }

    if (is.null(mvr_pars)) {
        cov_mat <- NULL
    }
    return(list(
        "chain_file" = mcmc_chain_file, "history_file" = infection_history_file,
        "cov_mat" = cov_mat, "step_scale" = steps,
        "overall_swap_proposals"=overall_swap_proposals,
        "overall_add_proposals"=overall_add_proposals
        
    ))
}
