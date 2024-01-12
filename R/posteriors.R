

#' Posterior function pointer
#'
#' Takes all of the input data/parameters and returns a function pointer. This function finds the posterior for a given set of input parameters (theta) and infection histories without needing to pass the data set back and forth. No example is provided for function_type=2, as this should only be called within \code{\link{serosolver}}
#' @param par_tab the parameter table controlling information such as bounds, initial values etc. See \code{\link{example_par_tab}}
#' @param antibody_data the data frame of data to be fitted. Must have columns: group (index of group); individual (integer ID of individual); samples (numeric time of sample taken); virus (numeric time of when the virus was circulating); biomarker_group (integer of the observation group type, using a unique value for each distinctive type of observation underpinned by the same generative model); titre (integer of titre value against the given virus at that sampling time). See \code{\link{example_antibody_data}}
#' @param antigenic_map (optional) a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. See \code{\link{example_antigenic_map}}
#' @param possible_exposure_times (optional) if no antigenic map is specified, this argument gives the vector of times at which individuals can be infected
#' @param prior_version which infection history assumption version to use? See \code{\link{describe_priors}} for options. Can be 1, 2, 3 or 4
#' @param solve_likelihood usually set to TRUE. If FALSE, does not solve the likelihood and instead just samples/solves based on the model prior
#' @param age_mask see \code{\link{create_age_mask}} - a vector with one entry for each individual specifying the first epoch of circulation in which an individual could have been exposed
#' @param measurement_indices_by_time if not NULL, then use these indices to specify which measurement bias parameter index corresponds to which time
#' @param n_alive if not NULL, uses this as the number alive in a given year rather than calculating from the ages. This is needed if the number of alive individuals is known, but individual birth dates are not
#' @param function_type integer specifying which version of this function to use. Specify 1 to give a posterior solving function; 2 to give the gibbs sampler for infection history proposals; otherwise just solves the titre model and returns predicted titres. NOTE that this is not the same as the attack rate prior argument, \code{version}!
#' @param titre_before_infection TRUE/FALSE value. If TRUE, solves titre predictions, but gives the predicted titre at a given time point BEFORE any infection during that time occurs.
#' @param data_type integer, currently accepting 1 for discrete or 2 for continuous. 
#' @param VERBOSE if TRUE, prints warning messages
#' @param ... other arguments to pass to the posterior solving function
#' @return a single function pointer that takes only pars and infection_histories as unnamed arguments. This function goes on to return a vector of posterior values for each individual
#' @examples
#' \dontrun{
#' data(example_par_tab)
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' data(example_inf_hist)
#'
#' ## Simple model solving code. Output matches entries of example_antibody_data
#' model_func <- create_posterior_func(example_par_tab, example_antibody_data, example_antigenic_map, function_type = 3)
#' y <- model_func(example_par_tab$values, example_inf_hist)
#'
#' ## Solve likelihood
#' par_tab <- example_par_tab[example_par_tab$names != "phi",]
#' likelihood_func <- create_posterior_func(par_tab, example_antibody_data, example_antigenic_map, function_type = 1, prior_version = 2)
#' liks <- likelihood_func(par_tab$values, example_inf_hist)
#' }
#' @export
create_posterior_func <- function(par_tab,
                                  antibody_data,
                                  antigenic_map=NULL,
                                  possible_exposure_times=NULL,
                                  prior_version = 1,
                                  solve_likelihood = TRUE,
                                  age_mask = NULL,
                                  measurement_indices_by_time = NULL,
                                  n_alive = NULL,
                                  function_type = 1,
                                  antibody_level_before_infection=FALSE,
                                  data_type=1,
                                  biomarker_groups_weights =1,
                                  VERBOSE=FALSE,
                                  ...) {
    check_par_tab(par_tab, TRUE, prior_version,VERBOSE)
    if (!("population_group" %in% colnames(antibody_data))) {
        antibody_data$population_group <- 1
    }
   
    ## Add a dummy observation type variable if not provided
    if (!("biomarker_group" %in% colnames(antibody_data))) {
        if(VERBOSE) message(cat("Note: no biomarker_group detected in antibody_data. Assuming all biomarker_group as 1."))
        antibody_data$biomarker_group <- 1
    }
    
    if (!("biomarker_group" %in% colnames(par_tab))) {
      if(VERBOSE) message(cat("Note: no biomarker_group detected in par_tab. Assuming all biomarker_group as 1."))
        par_tab$biomarker_group <- 1
    }
    if (!is.null(measurement_indices_by_time) & !("biomarker_group" %in% colnames(measurement_indices_by_time))) {
      if(VERBOSE) message(cat("Note: no biomarker_group detected in measurement_indices_by_time. Assuming all biomarker_group as 1."))
      measurement_indices_by_time$biomarker_group <- 1
    }
  
    ## Check that antibody data is formatted correctly
    check_data(antibody_data,VERBOSE)
    antibody_data <- antibody_data %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
    
    ## Get unique observation types
    unique_biomarker_groups <- unique(antibody_data$biomarker_group)
    unique_biomarker_groups <- unique_biomarker_groups[order(unique_biomarker_groups)]
    n_biomarker_groups <- length(unique_biomarker_groups)
    
    n_indivs <- length(unique(antibody_data$individual))
    
    ## Likelihood versions for different obs types
    if(length(data_type) ==1 & n_biomarker_groups > 1){
        data_type <- rep(data_type, n_biomarker_groups)
    }
    
    if(length(biomarker_groups_weights) ==1 & n_biomarker_groups > 1){
        biomarker_groups_weights <- rep(1, n_biomarker_groups)
    }
    
    #########################################################
    ## SETUP ANTIGENIC MAP
    #########################################################
    ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
    if (!is.null(antigenic_map)) {
        possible_exposure_times_tmp <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
        if(!is.null(possible_exposure_times) & !identical(possible_exposure_times, possible_exposure_times_tmp)){
          if(VERBOSE) message(cat("Warning: provided possible_exposure_times argument does not match entries in the antigenic map. Please make sure that there is an entry in the antigenic map for each possible circulation time. Using the antigenic map times."))
        }
      possible_exposure_times <- possible_exposure_times_tmp
      
      ## If no observation types assumed, set all to 1.
      if (!("biomarker_group" %in% colnames(antigenic_map))) {
        if(VERBOSE) message(cat("Note: no biomarker_group detection in antigenic_map. Aligning antigenic map with par_tab."))
          antigenic_map_tmp <- replicate(n_biomarker_groups,antigenic_map,simplify=FALSE)
          for(biomarker_group in unique_biomarker_groups){
              antigenic_map_tmp[[biomarker_group]]$biomarker_group <- biomarker_group
          }
          antigenic_map <- do.call(rbind,antigenic_map_tmp)
      }
      
    } else {
        ## Create a dummy map with entries for each observation type
      antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,
                                  "inf_times"=rep(possible_exposure_times, n_biomarker_groups), 
                                  "biomarker_group"=rep(unique_biomarker_groups,each=length(possible_exposure_times)))
    }
    #########################################################
    ## SETUP DATA
    #########################################################
    ## Separate out initial readings and repeat readings - we only
    ## want to solve the model once for each unique indiv/sample/biomarker_id year tested
    antibody_data_unique <- antibody_data[antibody_data$repeat_number == 1, ]
    ## Observations from repeats
    antibody_data_repeats <- antibody_data[antibody_data$repeat_number != 1, ]
    ## Find which entry in antibody_data_unique each antibody_data_repeats entry should correspond to
    tmp <- row.match(
        antibody_data_repeats[, c("individual", "sample_time", "biomarker_group", "biomarker_id")],
        antibody_data_unique[, c("individual", "sample_time", "biomarker_group", "biomarker_id")]
    )
    antibody_data_repeats$index <- tmp

    ## Which entries in the overall antibody_data matrix does each entry in antibody_data_unique correspond to?
    overall_indices <- row.match(
        antibody_data[, c("individual", "sample_time", "biomarker_group","biomarker_id")],
        antibody_data_unique[, c("individual", "sample_time", "biomarker_group","biomarker_id")]
    )

    ## Setup data vectors and extract
    setup_dat <- setup_antibody_data_for_posterior_func(
        antibody_data_unique, antigenic_map, 
        possible_exposure_times,
        age_mask, n_alive
    )
    ## Vector of observation types matching the unique samples
    biomarker_groups <- setup_dat$biomarker_groups
    
    ## Number of unique groups
    n_groups <- length(unique(antibody_data$population_group))
    group_id_vec <- setup_dat$group_id_vec
    
    ## List of melted antigenic maps, one entry for each observation type
    antigenic_map_melted <- setup_dat$antigenic_map_melted
    antigenic_distances <- antigenic_map_melted[[1]]
    
    possible_exposure_times <- setup_dat$possible_exposure_times
    exposure_id_indices <- setup_dat$exposure_id_indices
    
    ## Sample collection times, entry for each unique individual, observation type and sample
    sample_times <- setup_dat$sample_time
    ## Indices related to entries in sample_data
    sample_data_start <- setup_dat$sample_data_start
    
    ## Indices related to entries in antibody_data
    nrows_per_sample <- setup_dat$nrows_per_sample
    antibody_data_start <- setup_dat$antibody_data_start
    
    ## Indices related to entries in type_data
    type_data_start <- setup_dat$type_data_start
    biomarker_groups <- setup_dat$biomarker_groups
    
    ## Indices related to entries in the antigenic map
    biomarker_id_indices <- setup_dat$biomarker_id_indices
    
    n_alive <- setup_dat$n_alive
    age_mask <- setup_dat$age_mask
    sample_mask <- setup_dat$sample_mask
    n_indiv <- setup_dat$n_indiv
    DOBs <- setup_dat$DOBs

    ## Which entries of the unique antibody data correspond to each individual? 
    ## Used to summarize into per-individual likelihoods later
    nrows_per_individual_in_data <-  antibody_data_unique %>% group_by(individual, biomarker_group) %>% 
        tally() %>% 
        tidyr::pivot_wider(id_cols=individual,values_from=n,names_from=biomarker_group,values_fill=0) %>%
        ungroup() %>%
        select(-individual)
    nrows_per_individual_in_data <- nrows_per_individual_in_data %>% 
        select(order(colnames(nrows_per_individual_in_data)))
    nrows_per_individual_in_data <- as.matrix(nrows_per_individual_in_data)
    
    tmp <- antibody_data_unique %>% group_by(individual, biomarker_group) %>% 
        tally() %>% pull(n)
    cum_nrows_per_individual_in_data <- cumsum(c(0,tmp))
    
    ## Some additional setup for the repeat data
    ## Used to summarize into per-individual likelihoods later
    nrows_per_individual_in_data_repeats <-  antibody_data_repeats %>% group_by(individual, biomarker_group) %>% 
        tally() %>% 
        tidyr::pivot_wider(id_cols=individual,values_from=n,names_from=biomarker_group,values_fill=0) %>%
        ungroup()
    
    ## Need to make sure that there are entries for each individual, even if 0s
    indiv_repeat_indices <- nrows_per_individual_in_data_repeats %>% pull(individual)
    
    nrows_per_individual_in_data_repeats <- nrows_per_individual_in_data_repeats %>%
        select(-individual)
    nrows_per_individual_in_data_repeats <- nrows_per_individual_in_data_repeats %>% 
        select(order(colnames(nrows_per_individual_in_data_repeats)))
    nrows_per_individual_in_data_repeats <- as.matrix(nrows_per_individual_in_data_repeats)
    
    tmp <- antibody_data_repeats %>% group_by(individual, biomarker_group) %>% 
        tally() %>% pull(n)
    cum_nrows_per_individual_in_data_repeats <- cumsum(c(0,tmp))
    
    biomarker_group_indices <- lapply(unique_biomarker_groups, function(x) which(antibody_data_unique$biomarker_group == x))
    biomarker_group_indices_repeats <- lapply(unique_biomarker_groups, function(x) which(antibody_data_repeats$biomarker_group == x))
    
    ## Pull out unique and repeat antibody levels for solving likelihood later
    antibody_levels_unique <- antibody_data_unique$measurement
    antibody_levels_repeats <- antibody_data_repeats$measurement
    repeat_indices <- antibody_data_repeats$index
    repeat_indices_cpp <- repeat_indices - 1
    
    ## Not a good solution, but make lists of the antibody_level data, repeat antibody_level data and indices. This will be used by the Cpp proposal
    ## function
    antibody_levels_unique_list <- list()
    antibody_levels_repeats_list <- list()
    repeat_indices_cpp_list <- list()
    for(x in unique_biomarker_groups){
        tmp_antibody_data <- antibody_data_unique[antibody_data_unique$biomarker_group == x,]
        tmp_antibody_data_repeats <- antibody_data_repeats[antibody_data_repeats$biomarker_group == x,]
        
        tmp <- row.match(
            tmp_antibody_data_repeats[, c("individual", "sample_time", "biomarker_group", "biomarker_id")],
            tmp_antibody_data[, c("individual", "sample_time", "biomarker_group", "biomarker_id")]
        )
        
        antibody_levels_unique_list[[x]] <- tmp_antibody_data$measurement
        antibody_levels_repeats_list[[x]] <- tmp_antibody_data_repeats$measurement
        repeat_indices_cpp_list[[x]] <- tmp-1
    }
    
    #########################################################
    ## PARAMETER TABLE
    #########################################################
    ## Extract parameter type indices from par_tab, to split up
    ## similar parameters in model solving functions
    
    ## In general we are just going to use the indices for a single observation type
    par_tab_unique <- par_tab[!is.na(par_tab$biomarker_group) & par_tab$biomarker_group == min(par_tab$biomarker_group),]

    ## These will be different for each biomarker_group
    option_indices <- which(par_tab$par_type == 0)
    theta_indices <- which(par_tab$par_type %in% c(0, 1))
    measurement_indices_par_tab <- which(par_tab$par_type == 3)

    theta_indices_unique <- which(par_tab_unique$par_type %in% c(0, 1))
    
    ## Each biomarker_group must have the same vector of parameters in the same order
    par_names_theta <- par_tab_unique[theta_indices_unique, "names"]
    theta_indices_unique <- seq_along(theta_indices_unique) - 1
    names(theta_indices_unique) <- par_names_theta
    
    par_names_theta_all <- par_tab[theta_indices,"names"]
    n_pars <- length(theta_indices_unique)
    
    ## Sort out any assumptions for measurement bias
    use_measurement_bias <- (length(measurement_indices_par_tab) > 0) & !is.null(measurement_indices_by_time)
    
    antibody_level_shifts <- c(0)
    expected_indices <- NULL
    measurement_bias <- NULL
    additional_arguments <- NULL

    repeat_data_exist <- nrow(antibody_data_repeats) > 0
    if (use_measurement_bias) {
        if(VERBOSE) message(cat("Using measurement bias\n"))
        expected_indices <- antibody_data_unique %>% left_join(measurement_indices_by_time,by = c("biomarker_id", "biomarker_group")) %>% pull(rho_index)
    } else {
        expected_indices <- c(-1)
    }
  
    repeat_indices_bool <- TRUE
    if (!repeat_data_exist) {
        repeat_indices_cpp <- c(-1)
        repeat_indices_bool <- FALSE
    }

    ## These will be the same for each biomarker_group, as currently only one exposure type
    phi_indices <- which(par_tab$par_type == 2)
    ## weights_indices <- which(par_tab$par_type == 4) ## For functional form version of FOI
    ## knot_indices <- which(par_tab$par_type == 5) ## For functional form version of FOI
    
    ## Find which options are being used in advance for speed
    explicit_phi <- "phi" %in% par_tab$names ## Explicit prob of infection term
    # spline_phi <- (length(knot_indices) > 0)
    
    ## Allow different likelihood function for each observation type
    ## Set data type for likelihood function
    ## Potentially different likelihood for each observation type
    likelihood_func_use <- list()
    for(biomarker_group in unique_biomarker_groups){
        if(data_type[biomarker_group] == 1){
          likelihood_func_use[[biomarker_group]] <- likelihood_func_fast
          if(VERBOSE) message(cat("Setting to discretized, bounded observations\n"))
          
        } else if(data_type[biomarker_group] == 2){
          if(VERBOSE) message(cat("Setting to continuous, bounded observations\n"))
          likelihood_func_use[[biomarker_group]] <- likelihood_func_fast_continuous
        } else {
          if(VERBOSE) message(cat("Assuming discretized, bounded observations\n"))
          likelihood_func_use[[biomarker_group]] <- likelihood_func_fast
        }
    }
    
    if (function_type == 1) {
      if(VERBOSE) message(cat("Creating posterior solving function...\n"))
        f <- function(pars, infection_history_mat) {
          
          ## Transmission prob is the part of the likelihood function corresponding to each individual
          transmission_prob <- rep(0, n_indiv)
          if (explicit_phi) {
            phis <- pars[phi_indices]
            transmission_prob <- calc_phi_probs_indiv(
              phis, infection_history_mat,
              age_mask, sample_mask)
          }
          if (solve_likelihood) {
            theta <- pars[theta_indices]
            names(theta) <- par_names_theta_all
            
            antigenic_map_long <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
            antigenic_map_short <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
            
            cr_longs <- theta[which(par_names_theta_all=="cr_long")]
            cr_shorts <- theta[which(par_names_theta_all=="cr_short")]
            
            for(biomarker_group in unique_biomarker_groups){
              antigenic_map_long[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_longs[biomarker_group])
              antigenic_map_short[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_shorts[biomarker_group])
            }
            y_new <- antibody_model(
              theta, 
              theta_indices_unique, 
              unique_biomarker_groups,
              infection_history_mat, 
              possible_exposure_times, 
              exposure_id_indices,
              sample_times, 
              type_data_start,
              biomarker_groups,
              sample_data_start, 
              antibody_data_start,
              nrows_per_sample, 
              biomarker_id_indices, 
              antigenic_map_long,
              antigenic_map_short,
              antigenic_distances,
              antibody_level_before_infection
            )
            if (use_measurement_bias) {
              measurement_bias <- pars[measurement_indices_par_tab]
              antibody_level_shifts <- measurement_bias[expected_indices]
              y_new <- y_new + antibody_level_shifts
            }
                ## Calculate likelihood for unique antibody_levels and repeat data
                ## Sum these for each individual
                liks <- numeric(n_indivs)
                for(biomarker_group in unique_biomarker_groups){
                    ## Need theta for each observation type
                    liks_tmp <- likelihood_func_use[[biomarker_group]](
                                                    theta[(theta_indices_unique+1) + n_pars*(biomarker_group-1)], 
                                                    antibody_levels_unique[biomarker_group_indices[[biomarker_group]]], 
                                                    y_new[biomarker_group_indices[[biomarker_group]]])
                    
                    liks <- liks + biomarker_groups_weights[biomarker_group]*sum_buckets(liks_tmp, nrows_per_individual_in_data[,biomarker_group])
                    if (repeat_data_exist) {
                        ## Need theta for each observation type
                        
                        liks_repeats <- likelihood_func_use[[biomarker_group]](
                            theta[(theta_indices_unique+1) + n_pars*(biomarker_group-1)], 
                            antibody_levels_repeats[biomarker_group_indices_repeats[[biomarker_group]]], 
                            y_new[repeat_indices][biomarker_group_indices_repeats[[biomarker_group]]])
                        
                        liks[indiv_repeat_indices] <- liks[indiv_repeat_indices] + biomarker_groups_weights[biomarker_group]*sum_buckets(liks_repeats, nrows_per_individual_in_data_repeats[,biomarker_group])
                    }
                }
            } else {
                liks <- rep(-100000, n_indiv)
            }
            return(list(liks, transmission_prob))
        }
    } else if (function_type == 2) {
        
      if(VERBOSE) message(cat("Creating infection history proposal function\n"))
        if (prior_version == 4) {
            n_alive_total <- rowSums(n_alive)
        } else {
            n_alive_total <- c(-1, -1)
        }
        infection_model_prior_shape1 <- par_tab[par_tab$names == "infection_model_prior_shape1","values"]
        infection_model_prior_shape2 <- par_tab[par_tab$names == "infection_model_prior_shape2","values"]
        n_infected_group <- c(0, 0)
        ## Generate prior lookup table
        lookup_tab <- create_prior_lookup_groups(antibody_data, possible_exposure_times, infection_model_prior_shape1, infection_model_prior_shape2, n_alive)
        
        ## Use the original gibbs proposal function if no antibody_level immunity
        f <- function(pars, infection_history_mat,
                      probs, sampled_indivs,
                      infection_model_prior_shape1, 
                      infection_model_prior_shape2,
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
            names(theta) <- par_names_theta_all

            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_par_tab]
                antibody_level_shifts <- measurement_bias[expected_indices]
            }
            
            antigenic_map_long <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
            antigenic_map_short <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
            
            cr_longs <- theta[which(par_names_theta_all=="cr_long")]
            cr_shorts <- theta[which(par_names_theta_all=="cr_short")]
            
            for(biomarker_group in unique_biomarker_groups){
                antigenic_map_long[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_longs[biomarker_group])
                antigenic_map_short[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_shorts[biomarker_group])
            }

            n_infections <- sum_infections_by_group(infection_history_mat, group_id_vec, n_groups)
            if (prior_version == 4) n_infected_group <- rowSums(n_infections)
            ## Now pass to the C++ function
            res <- inf_hist_prop_prior_v2_and_v4(
                theta,
                theta_indices_unique, 
                unique_biomarker_groups,
                infection_history_mat,
                probs,
                sampled_indivs,
                n_infs,
                age_mask,
                sample_mask,
                n_alive,
                n_infections,
                n_infected_group,
                lookup_tab,
                swap_propn,
                swap_dist,
                propose_from_prior,
                infection_model_prior_shape1,
                infection_model_prior_shape2,
                possible_exposure_times,
                exposure_id_indices,
                sample_times,
                
                type_data_start,
                biomarker_groups,
                sample_data_start, 
                antibody_data_start,
                nrows_per_sample,
                
                cum_nrows_per_individual_in_data,
                cum_nrows_per_individual_in_data_repeats,

                group_id_vec,
                biomarker_id_indices,
                
                antigenic_map_long,
                antigenic_map_short,
                antigenic_distances,
                
                antibody_levels_unique,
                antibody_levels_repeats,
                
                length(antibody_levels_unique),
                
                repeat_indices_cpp,
                repeat_indices_bool,
                
                antibody_level_shifts,
                proposal_iter = proposal_iter,
                accepted_iter = accepted_iter,
                proposal_swap = proposal_swap,
                accepted_swap = accepted_swap,
                overall_swap_proposals,
                overall_add_proposals,
                proposal_ratios,
                n_alive_total,
                data_type,
                biomarker_groups_weights,
                temp,
                solve_likelihood
            )
            return(res)
        }
    } else {
      if(VERBOSE) message(cat("Creating model solving function...\n"))
        ## Final version is just the model solving function
        f <- function(pars, infection_history_mat) {
            theta <- pars[theta_indices]
            names(theta) <- par_names_theta_all

            antigenic_map_long <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
            antigenic_map_short <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
            
            cr_longs <- theta[which(par_names_theta_all=="cr_long")]
            cr_shorts <- theta[which(par_names_theta_all=="cr_short")]
            
            for(biomarker_group in unique_biomarker_groups){
                antigenic_map_long[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_longs[biomarker_group])
                antigenic_map_short[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_shorts[biomarker_group])
            }

            y_new <- antibody_model(
                theta, theta_indices_unique, unique_biomarker_groups,
                infection_history_mat, possible_exposure_times, 
                exposure_id_indices,
                sample_times, type_data_start,biomarker_groups,
                sample_data_start, antibody_data_start,
                nrows_per_sample, biomarker_id_indices, 
                antigenic_map_long,
                antigenic_map_short,
                antigenic_distances,
                antibody_level_before_infection
            )
            if (use_measurement_bias) {
                measurement_bias <- pars[measurement_indices_par_tab]
                antibody_level_shifts <- measurement_bias[expected_indices]
                y_new <- y_new + antibody_level_shifts
            }
            y_new[overall_indices]
        }
    }
    f
}
