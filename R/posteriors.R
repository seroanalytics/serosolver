

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
#' @param measurement_bias if not NULL, then use these indices to specify which measurement bias parameter index corresponds to which time
#' @param n_alive if not NULL, uses this as the number alive in a given year rather than calculating from the ages. This is needed if the number of alive individuals is known, but individual birth dates are not
#' @param function_type integer specifying which version of this function to use. Specify 1 to give a posterior solving function; 2 to give the gibbs sampler for infection history proposals; otherwise just solves the titre model and returns predicted titres. NOTE that this is not the same as the attack rate prior argument, \code{version}!
#' @param titre_before_infection TRUE/FALSE value. If TRUE, solves titre predictions, but gives the predicted titre at a given time point BEFORE any infection during that time occurs.
#' @param data_type integer or vector, with an entry for each unique data type in `antibody_data`. Set to 1 for discrete data (e.g., fold dilution) or 2 for continuous (e.g., ELISA optical density). 
#' @param biomarker_groups_weights integer or vector, giving a factor to multiply the log-likelihood contribution of this data type towards the overall likelihood.
#' @param start_level character, to tell the model how to treat initial antibody levels. This uses the observed data to either select starting values for each unique `individual`, `biomarker_id` and `biomarker_group` combination. See \code{\link{create_start_level_data}}. One of "min", "max", "mean", "median", or "full_random". Any other entry assumes all antibody starting levels are set to 0. Can also pass a tibble or data frame of starting levels matching the output of \code{\link{create_start_level_data}}.
#' @param start_level_randomize if TRUE, and data is discretized, then sets the starting antibody level to a random value between floor(x) and floor(x) + 1. Does nothing if using continuous data.
#' @param demographics if not NULL, then a tibble for each individual (1:n_indiv) giving demographic variable entries. Most importantly must include "birth" as the birth time. This is used if, for example, you have a stratification grouping in `par_tab`
#' @param verbose if TRUE, prints warning messages
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
                                  prior_version = 2,
                                  solve_likelihood = TRUE,
                                  age_mask = NULL,
                                  measurement_bias = NULL,
                                  n_alive = NULL,
                                  function_type = 1,
                                  antibody_level_before_infection=FALSE,
                                  data_type=1,
                                  biomarker_groups_weights =1,
                                  start_level = "other",
                                  start_level_randomize=FALSE,
                                  demographics=NULL,
                                  demographic_groups=NULL,
                                  fixed_inf_hists=NULL,
                                  verbose=FALSE,
                                  ...) {
    check_par_tab(par_tab, TRUE, prior_version,verbose)
    antibody_data <- as.data.frame(antibody_data)
    ## Add a dummy observation type variable if not provided
    if (!("biomarker_group" %in% colnames(antibody_data))) {
        if(verbose) message(cat("Note: no biomarker_group detected in antibody_data. Assuming all biomarker_group as 1.\n"))
        antibody_data$biomarker_group <- 1
    }
    
    if (!("biomarker_group" %in% colnames(par_tab))) {
      if(verbose) message(cat("Note: no biomarker_group detected in par_tab. Assuming all biomarker_group as 1.\n"))
        par_tab$biomarker_group <- 1
    }
    if (!is.null(measurement_bias) & !("biomarker_group" %in% colnames(measurement_bias))) {
      if(verbose) message(cat("Note: no biomarker_group detected in measurement_bias Assuming all biomarker_group as 1.\n"))
      measurement_bias$biomarker_group <- 1
    }
    if(verbose & any(class(start_level) == "character")){
      message(cat("Setting starting antibody levels based on data using command ", start_level, " and randomizing starting antibody levels set to ", start_level_randomize, "\n"))
    }
  
    if(verbose & !is.null(demographics)) message("Using time-varying demographic groupings for parameter stratification\n")
    
  
    ## Check that antibody data is formatted correctly
    check_data(antibody_data,verbose)
    antibody_data <- antibody_data %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
    
    ## Check demographics is formatted correctly
    if(!is.null(demographics)){
      check_demographics(demographics, par_tab, verbose)
    }
    
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
    antigenic_map_tmp <- setup_antigenic_map(antigenic_map, possible_exposure_times, n_biomarker_groups, unique_biomarker_groups,FALSE)
    antigenic_map <- antigenic_map_tmp$antigenic_map
    possible_exposure_times <- antigenic_map_tmp$possible_exposure_times
    infection_history_mat_indices <- antigenic_map_tmp$infection_history_mat_indices
   
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
    tmp <- get_demographic_groups(par_tab,antibody_data,demographics, demographic_groups)
    demographic_groups <- tmp$demographic_groups
    use_demographic_groups <- tmp$use_demographic_groups
    use_timevarying_demographics <- tmp$timevarying_demographics

    setup_dat <- setup_antibody_data_for_posterior_func(
      par_tab,
        antibody_data_unique, antigenic_map, 
        possible_exposure_times,
        age_mask, n_alive, verbose,use_demographic_groups,
        demographics
    )
    ## Vector of observation types matching the unique samples
    biomarker_groups <- setup_dat$biomarker_groups
    antibody_data_demo_index <- setup_dat$antibody_data_demo_group_index
    ## Number of unique groups
    demographics <- setup_dat$demographics

    indiv_pop_group_indices <- setup_dat$indiv_pop_group_indices
    unique_indiv_pop_group_indices1 <- unique(indiv_pop_group_indices)
    n_groups <- length(unique_indiv_pop_group_indices1[!is.na(unique_indiv_pop_group_indices1)])
    indiv_group_indices <- setup_dat$indiv_group_indices
    n_demographic_groups <- nrow(demographic_groups)
    demographics_groups <- setup_dat$demographics_groups
    
    ## List of melted antigenic maps, one entry for each observation type
    antigenic_map_melted <- setup_dat$antigenic_map_melted
    antigenic_distances <- antigenic_map_melted[[1]]
    
    possible_exposure_times <- setup_dat$possible_exposure_times
    exposure_id_indices <- setup_dat$exposure_id_indices
    possible_biomarker_ids <- setup_dat$possible_biomarker_ids
    infection_history_mat_indices <- setup_dat$infection_history_mat_indices
    
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
    
    ## Starting antibody levels
    births <- antibody_data_unique$birth
    ## If specifying starting levels based on a provided data frame or tibble, set them here. Otherwise, calculate them from the antibody data.
    if(class(start_level) %in% c("data.frame","tibble")){
      start_levels <- start_level 
      start_levels <- antibody_data %>% 
        left_join(start_levels %>% select(individual, biomarker_id, biomarker_group, starting_level, start_index) %>%distinct(),
                  by=c("individual","biomarker_group","biomarker_id"))
    } else {
      start_levels <- create_start_level_data(antibody_data,start_level,start_level_randomize) %>% 
        arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
    }
    start_levels <- start_levels %>% filter(repeat_number == 1)

    start_level_indices <- start_levels$start_index - 1
    start_antibody_levels <- start_levels %>% select(individual, biomarker_group,biomarker_id, starting_level) %>% distinct() %>% pull(starting_level)
    
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
    nrows_per_individual_in_data_repeats <-  antibody_data_repeats %>% 
      group_by(individual, biomarker_group) %>% 
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
    tmp <- antibody_data_repeats %>% 
      group_by(individual, biomarker_group) %>% 
        tally() %>%
      bind_rows(expand_grid(individual=unique(antibody_data$individual),
                            biomarker_group=unique(antibody_data$biomarker_group),
                            n=0)) %>%
      group_by(individual,biomarker_group) %>%
      dplyr::summarize(n=sum(n)) %>%
      pull(n)
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
    stratification_pars <- setup_stratification_table(par_tab, demographic_groups)
    scale_table <- stratification_pars[[1]]
    unique_groups <- 1:nrow(demographic_groups)
    ## Extract parameter type indices from par_tab, to split up
    ## similar parameters in model solving functions
    ## In general we are just going to use the indices for a single observation type
    par_tab_unique <- par_tab[!is.na(par_tab$biomarker_group) & par_tab$biomarker_group == min(par_tab$biomarker_group),]

    ## These will be different for each biomarker_group
    theta_indices <- which(par_tab$par_type %in% c(0, 1)) ## Which parameters are for the antibody kinetics model?

    scale_par_indices <- which(par_tab$par_type == 4) ## Which parameters are to scale the parameter values by groups?
    measurement_indices_par_tab <- which(par_tab$par_type == 3) ## Which parameters are measurement offsets?
    
    theta_meas_comb_indices <- c(theta_indices, measurement_indices_par_tab) #which(par_tab$par_type %in% c(0,1,3)) ## Both measurement offset and kinetics parameters
    
    theta_indices_unique <- which(par_tab_unique$par_type %in% c(0, 1)) ## Which parameters are for the antibody kinetics model, if we only had one biomarker group?
    rho_indices_unique <- which(par_tab_unique$par_type == 3)
    theta_meas_comb_indices_unique <- which(par_tab_unique$par_type %in% c(0,1,3))
    
    ## Find parameter transforms
    transforms <- par_tab[theta_meas_comb_indices,] %>% mutate(transform=case_when(
      lower_bound == 0 & upper_bound > 1 ~ 0, ## Log link
      lower_bound == 0 & upper_bound == 1 ~ 1, ## Logit link
      .default = 2  ## Otherwise no link
    )) %>% pull(transform)
    
    ## Each biomarker_group must have the same vector of parameters in the same order
    par_names_theta <- par_tab_unique[theta_indices_unique, "names"]
    par_names <- par_tab$names
    theta_indices_unique <- seq_along(theta_indices_unique) - 1
    names(theta_indices_unique) <- par_names_theta
    par_names_theta_all <- par_tab[theta_indices,"names"]
    n_pars <- length(theta_indices_unique)
    n_rhos <- length(measurement_indices_par_tab)
    
    ## Sort out any assumptions for measurement bias
    use_measurement_bias <- (length(measurement_indices_par_tab) > 0) & !is.null(measurement_bias)
    antibody_level_shifts <- c(0)
    expected_indices <- NULL
    measurement_bias_indices <- NULL
    additional_arguments <- NULL

    repeat_data_exist <- nrow(antibody_data_repeats) > 0
    if (use_measurement_bias) {
        if(verbose) message(cat("Using measurement bias\n"))
        ## Check that measurement bias is formatted correctly
        if(nrow(antibody_data_unique %>% select(biomarker_id,biomarker_group) %>% distinct()) > nrow(measurement_bias)) stop("not enough rho values provided in measurement_bias - one entry is needed for each unique combination of biomarker_id/biomarker_group")
      
        expected_indices <- antibody_data_unique %>% left_join(measurement_bias,by = c("biomarker_id", "biomarker_group")) %>% pull(rho_index)
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
          if(verbose) message(cat("Setting biomarker_group", biomarker_group, "to discretized, bounded observations\n"))
        } else if(data_type[biomarker_group] == 2){
          if(verbose) message(cat("Setting biomarker_group", biomarker_group, "to continuous, bounded observations\n"))
          likelihood_func_use[[biomarker_group]] <- likelihood_func_fast_continuous
        } else if(data_type[biomarker_group] == 3){
          if(verbose) message(cat("Setting biomarker_group", biomarker_group, "to continuous, bounded observations with false positives\n"))
          likelihood_func_use[[biomarker_group]] <- likelihood_func_fast_continuous_fp
        } else {
          if(verbose) message(cat("Assuming discretized, bounded observations for biomarker_group", biomarker_group, "\n"))
          likelihood_func_use[[biomarker_group]] <- likelihood_func_fast
        }
    }
    
    if (function_type == 1) {
      if(verbose) message(cat("Creating posterior solving function\n"))
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
            names(pars) <- par_names
            theta <- transform_parameters_cpp(pars, scale_table, theta_meas_comb_indices-1, scale_par_indices-1,as.matrix(demographic_groups), transforms)
            colnames(theta) <- names(pars[theta_meas_comb_indices])
            antigenic_map_long <- array(dim=c(length(possible_biomarker_ids)^2,n_biomarker_groups,n_demographic_groups))
            antigenic_map_short <- array(dim=c(length(possible_biomarker_ids)^2,n_biomarker_groups,n_demographic_groups))
            cr_longs <- matrix(theta[,colnames(theta)=="cr_long"],nrow=length(unique_groups))
            cr_shorts <- matrix(theta[,colnames(theta)=="cr_short"],nrow=length(unique_groups))
            #cr_shorts <- cr_shorts*cr_longs
            

            for(group in unique_groups){
              for(biomarker_group in unique_biomarker_groups){
                antigenic_map_long[,biomarker_group,group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_longs[group,biomarker_group])
                antigenic_map_short[,biomarker_group,group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_shorts[group,biomarker_group])
              }
            }
            y_new <- antibody_model(
              theta, 
              theta_indices_unique, 
              unique_biomarker_groups,
              infection_history_mat, 
              infection_history_mat_indices,
              indiv_group_indices,
              possible_exposure_times, 
              exposure_id_indices,
              sample_times, 
              type_data_start,
              biomarker_groups,
              sample_data_start, 
              antibody_data_start,
              nrows_per_sample, 
              biomarker_id_indices, 
              
              start_level_indices,
              start_antibody_levels,
              births,
              
              antigenic_map_long,
              antigenic_map_short,
              antigenic_distances,
              use_timevarying_demographics,
              antibody_level_before_infection
            )
            if (use_measurement_bias) {
              measurement_bias_indices <- matrix(theta[,colnames(theta) == "rho"],ncol=n_rhos)
              antibody_level_shifts <- measurement_bias_indices[cbind(antibody_data_demo_index,expected_indices)]
              y_new <- y_new + antibody_level_shifts
            }
                ## Calculate likelihood for unique antibody_levels and repeat data
                ## Sum these for each individual
                liks <- numeric(n_indivs)

                for(biomarker_group in unique_biomarker_groups){
                    ## Need theta for each observation type
                    liks_tmp <- likelihood_func_use[[biomarker_group]](
                      pars[theta_indices][(theta_indices_unique+1) + n_pars*(biomarker_group-1)], 
                                                    antibody_levels_unique[biomarker_group_indices[[biomarker_group]]], 
                                                    y_new[biomarker_group_indices[[biomarker_group]]])
                    
                    liks <- liks + biomarker_groups_weights[biomarker_group]*sum_buckets(liks_tmp, nrows_per_individual_in_data[,biomarker_group])
              
                    if (repeat_data_exist) {
                        ## Need theta for each observation type
                        
                        liks_repeats <- likelihood_func_use[[biomarker_group]](
                          pars[theta_indices][(theta_indices_unique+1) + n_pars*(biomarker_group-1)], 
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
      ## Enumerate out times each individual is able to be infected during
      indiv_possible_exposure_times <- tibble(
        individual = seq_along(age_mask),
        t_index = mapply(function(x, y) seq(x, y), age_mask, sample_mask, SIMPLIFY = FALSE)
      ) %>%
        unnest(cols = c(t_index)) %>%
        arrange(individual, t_index) %>%
        mutate(t_index = t_index - 1) %>%
        as.data.frame()
      ## Mask known infection states from estimation
      if(!is.null(fixed_inf_hists)) {
        if(verbose) message(cat("Fixing infection states given by fixed_inf_hists\n"))
        fixed_inf_hists$time <- match(fixed_inf_hists$time, possible_exposure_times)
        indiv_possible_exposure_times <- indiv_possible_exposure_times %>% 
          left_join(fixed_inf_hists %>% rename(t_index=time) %>% mutate(t_index =t_index - 1),
                    by=c("individual","t_index"))
        indiv_possible_exposure_times <- indiv_possible_exposure_times %>% filter(is.na(value))
      }
      
      ## Merge in any fixed infection states and remove these from consideration
      indiv_possible_exposure_times_indices <- indiv_possible_exposure_times %>% pull(t_index)
      ## Get start and end index in this vector for each individual
    indiv_poss_exp_times_start <- indiv_possible_exposure_times %>% 
      dplyr::mutate(i=1:n() - 1) %>%
      dplyr::group_by(individual) %>%
      dplyr::filter(t_index == min(t_index)) %>% 
      dplyr::ungroup() %>%
      complete(individual=1:n_indivs, fill=list(t_index=-1,i=-1,value=NA)) %>%
      dplyr::pull(i)
      
    indiv_poss_exp_times_end <- indiv_possible_exposure_times %>% 
      dplyr::mutate(i=1:n() - 1) %>%
      dplyr::group_by(individual) %>%
      dplyr::filter(t_index == max(t_index)) %>% 
      dplyr::ungroup() %>%
      complete(individual=1:n_indivs, fill=list(t_index=-1,i=-1,value=NA)) %>%
      dplyr::pull(i)
      
      if(verbose) message(cat("Creating infection history proposal function\n"))
        if (prior_version == 4) {
            n_alive_total <- rowSums(n_alive)
        } else {
            n_alive_total <- c(-1, -1)
        }
        infection_model_prior_shape1 <- par_tab[par_tab$names == "infection_model_prior_shape1","values"]
        infection_model_prior_shape2 <- par_tab[par_tab$names == "infection_model_prior_shape2","values"]
        n_infected_group <- rep(0, n_groups)
        ## Generate prior lookup table
        lookup_tab <- create_prior_lookup_groups(antibody_data, demographics,
                                                 possible_exposure_times[infection_history_mat_indices+1], 
                                                 infection_model_prior_shape1, infection_model_prior_shape2, n_alive)
        
        ## Use the original gibbs proposal function if no antibody_level immunity
        f <- function(pars, infection_history_mat,
                      probs, sampled_indivs,
                      infection_model_prior_shape1, 
                      infection_model_prior_shape2,
                      n_infs, proposal_inf_hist_indiv_swap_ratio,
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

          names(pars) <- par_names
          theta <- transform_parameters_cpp(pars, scale_table, theta_meas_comb_indices-1, scale_par_indices-1,as.matrix(demographic_groups), transforms)
          colnames(theta) <- names(pars[theta_meas_comb_indices])
          if (use_measurement_bias) {
            measurement_bias_indices <- matrix(theta[,colnames(theta) == "rho"],ncol=n_rhos)
            antibody_level_shifts <- measurement_bias_indices[cbind(antibody_data_demo_index,expected_indices)]
          }
          
          antigenic_map_long <- array(dim=c(length(possible_biomarker_ids)^2,n_biomarker_groups,n_demographic_groups))
          antigenic_map_short <- array(dim=c(length(possible_biomarker_ids)^2,n_biomarker_groups,n_demographic_groups))
          
          cr_longs <- matrix(theta[,colnames(theta)=="cr_long"],nrow=length(unique_groups))
          cr_shorts <- matrix(theta[,colnames(theta)=="cr_short"],nrow=length(unique_groups))
          #cr_shorts <- cr_shorts*cr_longs
          
          for(group in unique_groups){
            for(biomarker_group in unique_biomarker_groups){
              antigenic_map_long[,biomarker_group,group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_longs[group,biomarker_group])
              antigenic_map_short[,biomarker_group,group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_shorts[group,biomarker_group])
            }
          }
            n_infections <- sum_infections_by_group(infection_history_mat, indiv_pop_group_indices, n_groups,use_timevarying_demographics)
            n_infected_group <- rowSums(n_infections)
            #print("=========================================")
            #print("Theta: ")
            #print(theta)
            ## Now pass to the C++ function
            res <- inf_hist_prop_prior_v2_and_v4(
                theta,
                theta_indices_unique, 
                unique_biomarker_groups,
                indiv_group_indices,
                infection_history_mat,
                infection_history_mat_indices,
                probs,
                sampled_indivs,
                n_infs,
                age_mask,
                sample_mask,
                n_alive,
                n_infections,
                n_infected_group,
                lookup_tab,
                proposal_inf_hist_indiv_swap_ratio,
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

                indiv_pop_group_indices,
                biomarker_id_indices,
                
                start_level_indices,
                start_antibody_levels,
                births,
                
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
                indiv_possible_exposure_times_indices,
                indiv_poss_exp_times_start,
                indiv_poss_exp_times_end,
                use_timevarying_demographics,
                temp,
                solve_likelihood
            )
            return(res)
        }
    } else {
      if(verbose) message(cat("Creating model solving function\n"))
        ## Final version is just the model solving function
        f <- function(pars, infection_history_mat) {
          ## Need to create demographic-specific parameter transformations here
          names(pars) <- par_names
          theta <- transform_parameters_cpp(pars, scale_table, theta_meas_comb_indices-1, scale_par_indices-1,as.matrix(demographic_groups), transforms)
          colnames(theta) <- names(pars[theta_meas_comb_indices])
          antigenic_map_long <- array(dim=c(length(possible_biomarker_ids)^2,n_biomarker_groups,n_demographic_groups))
          antigenic_map_short <- array(dim=c(length(possible_biomarker_ids)^2,n_biomarker_groups,n_demographic_groups))
          
          cr_longs <- matrix(theta[,colnames(theta)=="cr_long"],nrow=length(unique_groups))
          cr_shorts <- matrix(theta[,colnames(theta)=="cr_short"],nrow=length(unique_groups))
          #cr_shorts <- cr_shorts*cr_longs
          for(group in unique_groups){
            for(biomarker_group in unique_biomarker_groups){
                antigenic_map_long[,biomarker_group,group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_longs[group,biomarker_group])
                antigenic_map_short[,biomarker_group,group] <- create_cross_reactivity_vector(antigenic_map_melted[[biomarker_group]], cr_shorts[group,biomarker_group])
            }
          }
          y_new <- antibody_model(
              theta, 
              theta_indices_unique, 
              unique_biomarker_groups,
              infection_history_mat, 
              infection_history_mat_indices,
              indiv_group_indices,
              
              possible_exposure_times, 
              exposure_id_indices,
              sample_times, 
              type_data_start,
              biomarker_groups,
              sample_data_start, 
              antibody_data_start,
              nrows_per_sample, 
              biomarker_id_indices, 
              
              start_level_indices,
              start_antibody_levels,
              births,
              
              antigenic_map_long,
              antigenic_map_short,
              antigenic_distances,
              use_timevarying_demographics,
              antibody_level_before_infection
          )

          if (use_measurement_bias) {
            measurement_bias_indices <- matrix(theta[,colnames(theta) == "rho"],ncol=n_rhos)
            antibody_level_shifts <- measurement_bias_indices[cbind(antibody_data_demo_index,expected_indices)]
            y_new <- y_new + antibody_level_shifts
          }
          y_new[overall_indices]
      }
    }
    f
}
