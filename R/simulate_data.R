#' Simulate full data set
#'
#' Simulates a full data set for a given set of parameters etc.
#' @param par_tab the full parameter table controlling parameter ranges and values
#' @param group which group index to give this simulated data
#' @param n_indiv number of individuals to simulate
#' @param antigenic_map (optional) A data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. See \code{\link{example_antigenic_map}}.
#' @param possible_exposure_times (optional) If no antigenic map is specified, this argument gives the vector of times at which individuals can be infected
#' @param measured_biomarker_ids vector of biomarker IDs that have titres measured matching entries in possible_exposure_times
#' @param sampling_times possible sampling times for the individuals, matching entries in possible_exposure_times
#' @param nsamps the number of samples each individual has (eg. nsamps=2 gives each individual 2 random sampling times from sampling_times)
#' @param missing_data numeric between 0 and 1, used to censor a proportion of titre observations at random (MAR)
#' @param age_min simulated age minimum
#' @param age_max simulated age maximum
#' @param attack_rates a vector of attack_rates for each entry in possible_exposure_times to be used in the simulation (between 0 and 1)
#' @param repeats number of repeat observations for each year
#' @param measurement_indices default NULL, optional vector giving the index of `measurement_bias` that each antigen/biomarker ID uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clusters, then this might be c(1,1,2,2,3,3)
#' @param obs_dist if not NULL, then a vector of data types to use for each biomarker_group
#' @return a list with: 1) the data frame of antibody data as returned by \code{\link{simulate_group}}; 2) a matrix of infection histories as returned by \code{\link{simulate_infection_histories}}; 3) a vector of ages
#' @family simulation_functions
#' @examples
#' data(example_par_tab)
#' data(example_antigenic_map)
#'
#' ## Times at which individuals can be infected
#' possible_exposure_times <- example_antigenic_map$inf_times
#' ## Simulate some random attack rates between 0 and 0.2
#' attack_rates <- runif(length(possible_exposure_times), 0, 0.2)
#' ## Vector giving the circulation times of measured antigens
#' sampled_antigens <- seq(min(possible_exposure_times), max(possible_exposure_times), by=2)
#' all_simulated_data <- simulate_data(par_tab=example_par_tab, group=1, n_indiv=50,    
#'                                    possible_exposure_times=possible_exposure_times,
#'                                    measured_biomarker_ids=sampled_antigens,
#'                                    sampling_times=2010:2015, nsamps=2, antigenic_map=example_antigenic_map, 
#'                                    age_min=10,age_max=75,
#'                                    attack_rates=attack_rates, repeats=2)
#' antibody_data <- all_simulated_data$data
#' antibody_data <- merge(antibody_data, all_simulated_data$ages)
#' @export
simulate_data <- function(par_tab,
                          group = 1,
                          n_indiv = 100,
                          antigenic_map = NULL,
                          possible_exposure_times = NULL,
                          measured_biomarker_ids = NULL,
                          sampling_times,
                          nsamps = 2,
                          missing_data = 0,
                          age_min = 5, age_max = 80,
                          attack_rates,
                          repeats = 1,
                          measurement_indices = NULL,
                          obs_dist = NULL) {

    #########################################################
    ## PARAMETER TABLE CHECKS
    #########################################################
    check_par_tab(par_tab)
    
    if (!("biomarker_group" %in% colnames(par_tab))) {
        message(cat("Note: no biomarker_group detection in par_tab Assuming all biomarker_group as 1. If this was deliberate, you can ignore this message."))
        par_tab$biomarker_group <- 1
    }
    
    ## Get unique observation types
    unique_biomarker_groups <- unique(par_tab$biomarker_group)
    n_biomarker_groups <- length(unique_biomarker_groups)
    
    #########################################################
    ## SETUP ANTIGENIC MAP
    #########################################################
    ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
    if (!is.null(antigenic_map)) {
        possible_exposure_times_tmp <- unique(antigenic_map$inf_times) # How many antigens are we testing against and what time did they circulate
        if(!is.null(possible_exposure_times) & !identical(possible_exposure_times, possible_exposure_times_tmp)){
            message(cat("Warning: provided possible_exposure_times argument does not match entries in the antigenic map. Please make sure that there is an entry in the antigenic map for each possible circulation time. Using the antigenic map times."))
        }
        possible_exposure_times <- possible_exposure_times_tmp
        
        ## If no observation types assumed, set all to 1.
        if (!("biomarker_group" %in% colnames(antigenic_map))) {
            message(cat("Note: no biomarker_group detection in antigenic_map. Aligning antigenic map with par_tab. If this was deliberate, you can ignore this message."))
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
    
    ## Check attack_rates entry
    check_attack_rates(attack_rates, possible_exposure_times)

    message("Simulating data\n")    

    
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

    ## Measurement indices unique
    measurement_indices_unique <- seq_along(which(par_tab_unique$par_type == 3)) - 1

    ## Extract parameters
    pars <- par_tab$values
    theta <- pars[theta_indices]
    names(theta) <- par_names_theta_all

    measurement_bias <- NULL
    if (!is.null(measurement_indices)) {
        message(cat("Measurement bias\n"))
        measurement_bias <- pars[measurement_indices_par_tab]
    }
    
    if (is.null(measured_biomarker_ids)) {
        measured_biomarker_ids <- possible_exposure_times
    }
    
    ## Simulate ages, where "age" is the age at the time of the last sample
    DOBs <- max(sampling_times) - floor(runif(n_indiv, age_min, age_max))
    
    ## Simulate infection histories
    tmp <- simulate_infection_histories(
        attack_rates, possible_exposure_times,
        sampling_times, DOBs
    )
    infection_history <- tmp[[1]]
    ARs <- tmp[[2]]
    ## Simulate antibody data
    sim_dat <- simulate_group(
        n_indiv, 
        
        theta, 
        theta_indices_unique,
        unique_biomarker_groups,
        
        infection_history,
        possible_exposure_times, 
        measured_biomarker_ids,
        sampling_times,
        nsamps, 
        antigenic_map, 
        repeats,
        measurement_bias,
        measurement_indices, 
        obs_dist, 
        DOBs
    )
    y <- sim_dat$antibody_data
    infection_history <- sim_dat$infection_history

    ## Need to update attack rate estimate based on strain mask, which is corrected in simulate_group
    age_mask <- create_age_mask(DOBs, possible_exposure_times)
    sample_mask <- create_sample_mask(y, possible_exposure_times)
    ## Can't be exposed or infected after the last sampling time
    for (i in 1:nrow(infection_history)) {
        if (sample_mask[i] < ncol(infection_history)) {
            infection_history[i, (sample_mask[i] + 1):ncol(infection_history)] <- 0
        }
    }

    n_alive <- sapply(seq(1, length(possible_exposure_times)), function(x)
        sum(age_mask <= x))
    ARs <- colSums(infection_history) / n_alive

    ## Randomly censor titre values
    y$measurement <- y$measurement * sample(c(NA, 1), nrow(y), prob = c(missing_data, 1 - missing_data), replace = TRUE)
    y$population_group <- group

    ages <- data.frame("individual" = 1:n_indiv, "birth" = DOBs)
    attack_rates <- data.frame("time" = possible_exposure_times, "AR" = ARs)
    return(list(
        data = y, infection_histories = infection_history,
        ages = ages, attack_rates = attack_rates, phis = attack_rates
    ))
}



#' Simulate group data
#'
#' Simulates a full set of antibody data for n_indiv individuals with known theta and infection_histories. Each individual gets nsamps random samples from sample_times, and infections can occur at any of possible_exposure_times
#' @inheritParams simulate_data
#' @param theta the named parameter vector
#' @param theta_indices_unique
#' @param unique_biomarker_groups
#' @param infection_histories the matrix of 1s and 0s giving presence/absence of infections for each individual
#' @param DOBs vector giving the time period of birth (entries corresponding to `possible_exposure_times`)
#' @return a data frame with columns individual, samples, virus and titre of simulated data
#' @family simulation_functions
#' @export
#' @seealso \code{\link{simulate_individual}}, \code{\link{simulate_individual_faster}}
simulate_group <- function(n_indiv,
                           theta,
                           theta_indices_unique,
                           unique_biomarker_groups,
                           infection_histories,
                           possible_exposure_times,
                           measured_biomarker_ids,
                           sample_times,
                           nsamps,
                           antigenic_map,
                           repeats = 1,
                           measurement_bias = NULL,
                           measurement_indices = NULL,
                           obs_dist = NULL,
                           DOBs=NULL) {
    n_biomarker_groups <- length(unique_biomarker_groups)
  ## Create antigenic map for short and long term boosting
    antigenic_map_long <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
    antigenic_map_short <- matrix(nrow=length(possible_exposure_times)^2, ncol=n_biomarker_groups)
    
    cr_longs <- theta[which(names(theta)=="cr_long")]
    cr_shorts <- theta[which(names(theta)=="cr_short")]
    
    for(biomarker_group in unique_biomarker_groups){
        antigenic_map_melted <- melt_antigenic_coords(antigenic_map[antigenic_map$biomarker_group == biomarker_group, c("x_coord", "y_coord")])
        antigenic_map_long[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted, cr_longs[biomarker_group])
        antigenic_map_short[,biomarker_group] <- create_cross_reactivity_vector(antigenic_map_melted, cr_shorts[biomarker_group])
    }
  antigenic_distances <- c(antigenic_map_melted)
  dat <- NULL
  ## For each individual
  for (i in 1:n_indiv) {
      if(!is.null(DOBs)) DOB <- DOBs[i]
    ## Choose random sampling times
    ## If there is one sampling time, then repeat the same sampling time
      sample_times_tmp <- sample_times[sample_times >= DOB]
      ## If all sample times are before birth, set to DOB
      if(length(sample_times_tmp) == 0){
          sample_times_tmp <- DOB
      }
    if (length(sample_times_tmp) == 1) {
      samps <- rep(sample_times_tmp, nsamps)
    } else {
      samps <- sample(sample_times_tmp, nsamps)
      samps <- samps[order(samps)]
    }

    ## Individuals can't be infected after their latest sampling time
    sample_mask <- max(which(max(samps) >= possible_exposure_times))
    if (sample_mask < ncol(infection_histories)) {
      infection_histories[i, (sample_mask + 1):ncol(infection_histories)] <- 0
    }
    y <- as.data.frame(simulate_individual_faster(
      theta,
      theta_indices_unique,
      infection_histories[i, ],
      antigenic_map_long,
      antigenic_map_short,
      antigenic_distances,
      unique_biomarker_groups,
      samps,
      possible_exposure_times,
      measured_biomarker_ids,
      measurement_bias,
      measurement_indices,
      obs_dist, repeats, DOB
    ))
    ## Record individual ID
    y$indiv <- i
    colnames(y) <- c("sample_time", "biomarker_id", "biomarker_group","measurement", "individual")
    ## Combine data
    dat <- rbind(dat, y[, c("individual", "sample_time", "biomarker_id", "biomarker_group","measurement")])
  }
  dat <- dat %>% group_by(individual,sample_time,biomarker_id,biomarker_group) %>% mutate(run = 1:n()) %>% ungroup() %>% as.data.frame()
  return(list(antibody_data = dat, infection_history = infection_histories))
}
#' Simulate individual data quickly
#'
#' FOR USERS: USE \code{\link{simulate_individual}}. This function does the same thing, but with a few short cuts for speed. Simulates a full set of antibody data for an individual with known theta and infection_history.
#' @inheritParams simulate_group
#' @param infection_history the vector of 1s and 0s giving presence/absence of infections
#' @param antigenic_map_long the long term antigenic cross reactivity map generated from \code{\link{create_cross_reactivity_vector}}
#' @param antigenic_map_short the short term antigenic cross reactivity map generated from \code{\link{create_cross_reactivity_vector}}
#' @param antigenic_distances (optional) same dimensions as antigenic_map_long and antigenic_map_short, but gives the raw euclidean antigenic distances
#' @param sampling_times vector of times at which blood samples were taken
#' @param measured_biomarker_ids vector of which biomarker IDs had measurements in `possible_exposure_times`
#' @return a data frame with columns samples, virus and titre of simulated data
#' @family simulation_functions
#' @export
simulate_individual_faster <- function(theta,
                                       unique_theta_indices=NULL,
                                       infection_history,
                                       antigenic_map_long,
                                       antigenic_map_short,
                                       antigenic_distances=NULL,
                                       unique_biomarker_groups=seq_len(ncol(antigenic_map_long)),
                                       sampling_times,
                                       possible_exposure_times,
                                       measured_biomarker_ids,
                                       measurement_bias = NULL, measurement_indices = NULL,
                                       obs_dist = NULL, 
                                       repeats = 1,
                                       DOB = NULL) {
 
    ## If have not passed the theta index vector, then create one based on unique parameter names
    if(is.null(unique_theta_indices)){
        unique_theta_names <- unique(names(thetas))
        unique_theta_indices <- seq_along(unique_theta_names)-1
        names(unique_theta_indices) <- unique_theta_names
    }
    n_pars <- length(unique_theta_indices)
  inf_hist <- matrix(nrow = 1, ncol = length(infection_history))
  inf_hist[1, ] <- infection_history

  n_samps <- length(sampling_times)
  n_biomarker_groups <- length(unique_biomarker_groups)
  
  sampling_times_long <- rep(sampling_times, n_biomarker_groups)
    
  ## length(measured_biomarker_ids) observations made per blood sample
  nrows_per_sample <- rep(length(measured_biomarker_ids), n_samps*n_biomarker_groups)

  ## Cumulative of the above for the algorithm
  antibody_data_start <- cumsum(c(0, nrows_per_sample))

  ## Iterate through sample times sample_times[0:(n_samps-1)] to solve the model
  sample_data_start <- cumsum(c(0, rep(n_samps,n_biomarker_groups)))

  ## Entries in the antigenic map
  strain_indices <- match(possible_exposure_times, possible_exposure_times) - 1

  ## Observation types to match sample times
  type_data_start <- c(0,n_biomarker_groups)
  
  ## Entries in the antigenic map for each measured strain
  measured_strain_indices <- match(rep(rep(measured_biomarker_ids, n_samps), n_biomarker_groups), possible_exposure_times) - 1
  dat <- matrix(nrow = length(measured_strain_indices) * repeats, ncol = 4) ## To store simulated data
  ## Go into C++ code to solve antibody model
  antibody_levels <- antibody_model(
    theta, 
    unique_theta_indices,
    unique_biomarker_groups,
    inf_hist, 
    possible_exposure_times, 
    strain_indices,
    sampling_times_long, 
    type_data_start=type_data_start,
    biomarker_groups=unique_biomarker_groups,
    sample_data_start, 
    antibody_data_start,
    nrows_per_sample, 
    measured_strain_indices,
    antigenic_map_long,
    antigenic_map_short,
    antigenic_distances,
    FALSE
  )
  ## Repeated each simulated titre per observation repeat
  antibody_levels <- rep(antibody_levels, repeats)
  ## Housekeeping for return data
  sampling_times <- rep(rep(sampling_times, each=length(measured_biomarker_ids)),n_biomarker_groups)
  #enum_repeats <- rep(1:repeats, each = length(sampling_times))
  sampling_times <- rep(sampling_times, repeats)
  dat[, 1] <- sampling_times
  dat[, 2] <- rep(rep(rep(measured_biomarker_ids, n_samps), repeats),n_biomarker_groups)
  
  biomarker_groups_data <- rep(unique_biomarker_groups,each=length(measured_biomarker_ids)*n_samps*repeats)
  
  dat[, 3] <- biomarker_groups_data
  ## Add observation noise, including measurement bias if selected
  if (!is.null(obs_dist)) {
      for(biomarker_group in unique_biomarker_groups){
        if (!is.null(measurement_indices)) {
            indices_tmp <- measurement_indices$biomarker_group == biomarker_group
            use_measurement_indices <- measurement_indices[indices_tmp,][match(dat[biomarker_groups_data==biomarker_group,2], measurement_indices[indices_tmp,"biomarker_id"]),"rho_index"]
          dat[biomarker_groups_data == biomarker_group, 4] <- add_noise(antibody_levels[biomarker_groups_data == biomarker_group], 
                                                                          theta[(unique_theta_indices+1) + n_pars*(biomarker_group-1)], 
                                                                          measurement_bias,use_measurement_indices,
                                                          data_type=obs_dist[biomarker_group])
        } else {
          dat[biomarker_groups_data == biomarker_group, 4] <- add_noise(antibody_levels[biomarker_groups_data == biomarker_group], 
                                                                          theta[(unique_theta_indices+1) + n_pars*(biomarker_group-1)], 
                                                                          NULL, NULL,data_type=obs_dist[biomarker_group])
        }
      }
  } else {
    dat[, 4] <- antibody_levels
  }
  return(dat)
}


#' Simulate individual data
#'
#' Simulates a full set of antibody data for an individual with known theta and infection_history.
#' @inheritParams simulate_group
#' @param infection_history the vector of 1s and 0s giving presence/absence of infections
#' @param sampling_times vector of times at which blood samples were taken
#' @param measured_biomarker_ids vector of which biomarker IDs had measurements in `possible_exposure_times`
#' @return a data frame with columns samples, biomarker IDs and antibody levels of simulated data
#' @family simulation_functions
#' @export
#' @examples
#' data(example_par_tab)
#' data(example_antigenic_map)
#' infection_history <- sample(c(0,1),size=nrow(example_antigenic_map), replace=TRUE,prob=c(0.9,0.1))
#' pars <- example_par_tab$values
#' names(pars) <- example_par_tab$names
#' possible_exposure_times <- example_antigenic_map$inf_times
#' y <- simulate_individual(pars, infection_history, example_antigenic_map, 2009, 
#'                          possible_exposure_times,possible_exposure_times,add_noise=FALSE)
simulate_individual <- function(theta,
                                infection_history,
                                antigenic_map,
                                sampling_times,
                                possible_exposure_times,
                                measured_biomarker_ids,
                                measurement_bias = NULL, measurement_indices = NULL,
                                add_noise = TRUE, repeats = 1,
                                DOB = NULL) {

  ## Create antigenic map for short and long term boosting
  antigenic_map_melted <- melt_antigenic_coords(antigenic_map[, c("x_coord", "y_coord")])
  antigenic_map_long <- create_cross_reactivity_vector(antigenic_map_melted, theta["cr_long"])
  antigenic_map_short <- create_cross_reactivity_vector(antigenic_map_melted, theta["cr_short"])
  antigenic_distances <- c(antigenic_map_melted)
  inf_hist <- matrix(nrow = 1, ncol = length(infection_history))
  inf_hist[1, ] <- infection_history

  n_samps <- length(sampling_times)
  ## length(measured_biomarker_ids) observatios made per blood sample
  rows_per_blood <- rep(length(measured_biomarker_ids), n_samps)

  ## Cumulative of the above for the algorithm
  cumu_rows <- c(0, sum(rows_per_blood))

  ## Iterate through sample times sample_times[0:(n_samps-1)] to solve the model
  rows_per_indiv <- c(0, n_samps)

  ## Entries in the antigenic map
  strain_indices <- match(possible_exposure_times, possible_exposure_times) - 1

  ## Entries in the antigenic map for each measured strain
  measured_biomarker_id_indices <- match(rep(measured_biomarker_ids, n_samps), possible_exposure_times) - 1
  dat <- matrix(nrow = length(measured_biomarker_ids) * repeats, ncol = 4) ## To store simulated data

  ## Go into C++ code to solve the antibody model
  antibody_levels <- antibody_model(
    theta, inf_hist, possible_exposure_times, strain_indices,
    sampling_times, rows_per_indiv, cumu_rows,
    rows_per_blood, measured_biomarker_ids,
    antigenic_map_long, antigenic_map_short,
    antigenic_distances
  )

  ## Repeated each simulated titre per observation repeat
  antibody_levels <- rep(antibody_levels, repeats)
  ## Housekeeping for return data
  sampling_times <- rep(sampling_times, rows_per_blood)
  enum_repeats <- rep(1:repeats, each = length(sampling_times))
  sampling_times <- rep(sampling_times, repeats)
  dat[, 1] <- sampling_times
  dat[, 2] <- rep(rep(measured_biomarker_ids, n_samps), repeats)
  ## Add observation noise, including measurement bias if selected
  if (add_noise) {
    if (!is.null(measurement_indices)) {
      dat[, 3] <- add_noise(antibody_levels, theta, measurement_bias, measurement_indices[match(dat[, 2], measured_biomarker_ids)])
    } else {
      dat[, 3] <- add_noise(antibody_levels, theta, NULL, NULL)
    }
  } else {
    dat[, 3] <- antibody_levels
  }
  dat[, 4] <- enum_repeats
  return(dat)
}

#' Add noise
#'
#' Adds truncated noise to antibody data
#' @param y the titre
#' @param theta a vector with max_measurement and error parameters
#' @param data_type integer, currently accepting 1 or 2. Set to 1 for discretized, bounded data, or 2 for continuous, bounded data. Note that with 2, min_measurement must be set.
#' @return a noisy titre
#' @export
#' @examples
#' \dontrun{
#' ## ... example in simulate_individual
#' pars <- c("obs_sd"=1)
#' y <- runif(100)
#' noisy_y <- add_noise(y, pars)
#' }
add_noise <- function(y, theta, measurement_bias = NULL, indices = NULL,data_type=1) {
  if(data_type ==2){
    if (!is.null(measurement_bias)) {
      noise_y <- rnorm(length(y), mean = y + measurement_bias[indices], sd = theta["obs_sd"])
    } else {
      noise_y <- rnorm(length(y), mean = y, sd = theta["obs_sd"])
    }
    
    ## If outside of bounds, truncate
    noise_y[noise_y < theta["min_measurement"]] <- theta["min_measurement"]
    noise_y[noise_y > theta["max_measurement"]] <- theta["max_measurement"]
  } else {
  ## Draw from normal
    if (!is.null(measurement_bias)) {
      noise_y <- floor(rnorm(length(y), mean = y + measurement_bias[indices], sd = theta["obs_sd"]))
    } else {
      noise_y <- floor(rnorm(length(y), mean = y, sd = theta["obs_sd"]))
    }
  
    ## If outside of bounds, truncate
    noise_y[noise_y < 0] <- 0
    noise_y[noise_y > theta["max_measurement"]] <- theta["max_measurement"]
  }
  return(noise_y)
}

#' Simulate attack rates
#'
#' Given a number of possible infection years, simulates attack rates from a log normal distribution with specified mean and standard deviation.
#' @param infection_years the number of infection years
#' @param mean_par the mean of the log normal
#' @param sd_par the sd of the log normal
#' @param large_first_year simulate an extra large attach rate in the first year?
#' @param big_year_mean if large first year, what mean to use?
#' @return a vector of attack rates
#' @family simulation_functions
#' @export
simulate_attack_rates <- function(infection_years, mean_par = 0.15, sd_par = 0.5,
                                  large_first_year = FALSE, big_year_mean = 0.5) {
  attack_year <- rlnorm(infection_years, meanlog = log(mean_par) - sd_par^2 / 2, sdlog = sd_par)
  if (large_first_year) attack_year[1] <- rlnorm(1, meanlog = log(big_year_mean) - (sd_par / 2)^2 / 2, sdlog = sd_par / 2)
  return(attack_year)
}

#' Simulate infection histories
#'
#' Given a vector of infection probabilities and potential infection times, simulates infections for each element of ages (ie. each element is an individual age. Only adds infections for alive individuals)
#' @param p_inf a vector of attack rates (infection probabilities) for each year
#' @param possible_exposure_times the vector of possible infection times
#' @param sampling_times vector of potential sampling times
#' @param DOBs a vector of ages for each individual
#' @return a list with a matrix of infection histories for each individual in ages and the true attack rate for each epoch
#' @family simulation_functions
#' @examples
#' p_inf <- runif(40,0.1,0.4)
#' possible_exposure_times <- seq_len(40) + 1967
#' n_indivs <- 100
#' sampling_times <- rep(max(possible_exposure_times), n_indivs)
#' DOBs <- rep(min(possible_exposure_times), n_indivs)
#' inf_hist <- simulate_infection_histories(p_inf, possible_exposure_times, sampling_times, DOBs)
#' @export
simulate_infection_histories <- function(p_inf, possible_exposure_times, sampling_times, DOBs) {
  n_strains <- length(p_inf) # How many strains
  n_indiv <- length(DOBs) # How many individuals
  indivs <- 1:n_indiv
  ## Empty matrix
  infection_histories <- matrix(0, ncol = n_strains, nrow = n_indiv)

  ## Simulate attack rates
  attack_rates <- p_inf

  ## Should this be necessary?
  attack_rates[attack_rates > 1] <- 1
  ARs <- numeric(n_strains)

  age_mask <- create_age_mask(DOBs, possible_exposure_times)

  ## For each strain (ie. each infection year)
  for (i in 1:n_strains) {
    ## If there are strains circulating beyond the max sampling times, then alive==0
    if (max(sampling_times) >= possible_exposure_times[i]) {
      ## Find who was alive (all we need sampling_times for is its max value)
      alive <- which(age_mask <= i)

      ## Sample a number of infections for the alive individuals, and set these entries to 1
      y <- round(length(indivs[alive]) * attack_rates[i])
      # y <- rbinom(1, length(indivs[alive]),attack_rates[i])
      ARs[i] <- y / length(indivs[alive])
      x <- sample(indivs[alive], y)
      infection_histories[x, i] <- 1
    } else {
      ARs[i] <- 0
    }
  }
  return(list(infection_histories, ARs))
}


#' Generates attack rates from an SIR model with fixed beta/gamma, specified final attack rate and the number of time "buckets" to solve over ie. buckets=12 returns attack rates for 12 time periods
generate_ar_annual <- function(AR, buckets) {
  SIR_odes <- function(t, x, params) {
    S <- x[1]
    I <- x[2]
    R <- x[3]
    inc <- x[4]

    beta <- params[1]
    gamma <- params[2]
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    dinc <- beta * S * I
    list(c(dS, dI, dR, dinc))
  }
  R0 <- 1.2
  gamma <- 1 / 5
  beta <- R0 * gamma
  t <- seq(0, 360, by = 0.1)
  results <- as.data.frame(deSolve::ode(
    y = c(S = 1, I = 0.0001, R = 0, inc = 0),
    times = t, func = SIR_odes,
    parms = c(beta, gamma)
  ))
  incidence <- diff(results$inc)
  incidence <- incidence * AR / sum(incidence)
  group <- 360 * 10 / buckets
  monthly_risk <- colSums(matrix(incidence, nrow = group))
  return(monthly_risk)
}


simulate_ars_buckets <- function(infection_years, buckets, mean_par = 0.15, sd_par = 0.5,
                                 large_first_year = FALSE, big_year_mean = 0.5) {
  n <- ceiling(length(infection_years) / buckets)
  attack_year <- rlnorm(n, meanlog = log(mean_par) - sd_par^2 / 2, sdlog = sd_par)
  if (large_first_year) attack_year[1] <- rlnorm(1, meanlog = log(big_year_mean) - (sd_par / 2)^2 / 2, sdlog = sd_par / 2)
  ars <- NULL

  for (i in seq_along(attack_year)) {
    ars <- c(ars, generate_ar_annual(attack_year[i], buckets))
  }

  ars <- ars[1:length(infection_years)]
  return(ars)
}

simulate_ars_spline <- function(infection_years, buckets, mean_par = 0.15, sd_par = 0.5, large_first_year = FALSE, big_year_mean = 0.5, knots, theta) {
  infection_years <- infection_years[seq(1, length(infection_years), by = buckets)] / buckets
  n <- length(infection_years)
  attack_year <- rlnorm(n, meanlog = log(mean_par) - sd_par^2 / 2, sdlog = sd_par)
  if (large_first_year) attack_year[1] <- rlnorm(1, meanlog = log(big_year_mean) - (sd_par / 2)^2 / 2, sdlog = sd_par / 2)
  ars <- generate_phis(attack_year, knots, theta, n, buckets)
  return(ars)
}
