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
#' @param data_type if not NULL, then a vector of data types to use for each biomarker_group
#' @param demographics if not NULL, then a tibble for each individual (1:n_indiv) giving demographic variable entries. Most importantly must include "birth" as the birth time. This is used if, for example, you have a stratification grouping in `par_tab`
#' @param verbose if TRUE, prints additional messages
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
                          age_group_bounds = NULL,
                          attack_rates,
                          repeats = 1,
                          measurement_indices = NULL,
                          data_type = NULL,
                          demographics=NULL,
                          verbose=FALSE) {
    #########################################################
    ## CHECK FOR BIOMARKER GROUPS
    #########################################################
    if (!("biomarker_group" %in% colnames(par_tab))) {
        if(verbose) message(cat("Note: no biomarker_group detection in par_tab Assuming all biomarker_group as 1. If this was deliberate, you can ignore this message.\n"))
        par_tab$biomarker_group <- 1
    }

    ## Get unique observation types
    unique_biomarker_groups <- unique(par_tab$biomarker_group)
    n_biomarker_groups <- length(unique_biomarker_groups)
    
    #########################################################
    ## SETUP ANTIGENIC MAP
    #########################################################
    antigenic_map_tmp <- setup_antigenic_map(antigenic_map, possible_exposure_times, 
                                             n_biomarker_groups,unique_biomarker_groups,FALSE)
    antigenic_map <- antigenic_map_tmp$antigenic_map
    possible_exposure_times <- antigenic_map_tmp$possible_exposure_times
    infection_history_mat_indices <- antigenic_map_tmp$infection_history_mat_indices
    ## Check attack_rates entry
    check_attack_rates(attack_rates, possible_exposure_times)
    
    ## Which biomarkers were measured
    if (is.null(measured_biomarker_ids)) {
      measured_biomarker_ids <- possible_exposure_times
    }

    #########################################################
    ## BUILD ANTIBODY DATA OBJECT
    #########################################################
    ## Create sampling times and birth dates
    if(is.null(demographics)){
      birth_dates <- tibble(individual=1:n_indiv,
                            birth=max(sampling_times) - floor(runif(n_indiv, age_min,age_max + 1)))
    } else {
       birth_dates <- demographics
    }
    sampling_times_data <- expand_grid(individual=1:n_indiv,sample_time=sampling_times) %>%
      left_join(birth_dates,by="individual") %>%
      ## Ensure only sampling after birth
      filter(sample_time >= birth) %>% 
      group_by(individual) %>%
      ## Ensure not sampling more than there are time periods to sample from
      sample_n(pmin(n(),nsamps))
    
    ## Expand out with measured biomarkers and repeats
    antibody_data <- expand_grid(sampling_times_data,repeat_number=1:repeats,
                                 biomarker_id = measured_biomarker_ids,
                                 biomarker_group=unique_biomarker_groups,
                                 measurement=0)
    
    antibody_data <- as.data.frame(antibody_data)
    
    ## If we've specified age_group_bounds, then we want to have timevarying age groups
    timevarying_demographics <- NULL
    if(!is.null(age_group_bounds)){
      timevarying_demographics <- demographics %>% expand_grid(time=possible_exposure_times)
      timevarying_demographics$age <- timevarying_demographics$time - timevarying_demographics$birth
      timevarying_demographics$age_group <- as.numeric(cut(timevarying_demographics$age, breaks=c(0, age_group_bounds)))
    }
    #########################################################
    ## PARAMETER TABLE CHECKS
    #########################################################
    if(!(4 %in% unique(par_tab$par_type))){
      par_tab <- add_scale_pars(par_tab,antibody_data, timevarying_demographics)
    }
    par_tab <- check_par_tab(par_tab)
    #########################################################
    ## SIMULATE DATA
    #########################################################
    message("Simulating data\n")    
    
    measurement_bias <- NULL
    if (!is.null(measurement_indices)) {
        message(cat("Measurement bias\n"))
        #measurement_bias <- pars[measurement_indices_par_tab]
    }
    ## Simulate infection histories
    DOBs <- antibody_data %>% dplyr::select(individual,birth) %>% distinct() %>% pull(birth)
    tmp <- simulate_infection_histories(
        attack_rates, possible_exposure_times,
        sampling_times, DOBs
    )
    infection_history <- tmp[[1]]
    ARs <- tmp[[2]]
   
    ## Need to update attack rate estimate based on strain mask, which is corrected in simulate_group
    age_mask <- create_age_mask(DOBs, possible_exposure_times)
    sample_mask <- create_sample_mask(antibody_data, possible_exposure_times)
    ## Can't be exposed or infected after the last sampling time
    for (i in 1:nrow(infection_history)) {
      if (sample_mask[i] < ncol(infection_history)) {
        infection_history[i, (sample_mask[i] + 1):ncol(infection_history)] <- 0
      }
    }
    
    ## Create model solving function
    ## Check that antibody data is formatted correctly
    check_data(antibody_data,verbose)
    antibody_data <- antibody_data %>% 
      arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
    f <- create_posterior_func(par_tab,antibody_data,antigenic_map,function_type=3,
                               possible_exposure_times = possible_exposure_times,
                               demographics=timevarying_demographics,
                               start_level="none")
    antibody_data$measurement <- f(par_tab$values, infection_history)
   
    ## Add noise, but need to be specific to the data type
    for(i in seq_along(unique_biomarker_groups)){
      ## Find indices for each unique biomarker group
      tmp_indices <- which(antibody_data$biomarker_group == unique_biomarker_groups[i])
      ## Get model parameters for this biomarker group
      tmp_par_tab <- par_tab[par_tab$biomarker_group == unique_biomarker_groups[i],]
      tmp_pars <- tmp_par_tab$values
      names(tmp_pars) <- tmp_par_tab$names
      ## Add noise, specific to the biomarker group's data type
      antibody_data$measurement[tmp_indices] <- add_noise(antibody_data$measurement[tmp_indices],tmp_pars,NULL,NULL,data_type=data_type[i])
     }
    ## Randomly censor titre values
    antibody_data <- antibody_data %>% mutate(measurement=if_else(runif(n())<missing_data,NA,measurement))
    antibody_data$population_group <- group
  n_alive <- get_n_alive(antibody_data,times = possible_exposure_times)
    ARs <- colSums(infection_history) / n_alive
    attack_rates_obs <- data.frame("time" = possible_exposure_times, "AR" = ARs)
    return(list(
        antibody_data = antibody_data, infection_histories = infection_history,
        attack_rates = attack_rates_obs, phis = attack_rates,par_tab=par_tab
    ))
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
    noise_y[noise_y < theta["min_measurement"]] <- theta["min_measurement"]
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
#' n_indiv <- 100
#' sampling_times <- rep(max(possible_exposure_times), n_indiv)
#' DOBs <- rep(min(possible_exposure_times), n_indiv)
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

#' Simulate the antibody model
#' 
#' Simulates the trajectory of the serosolver antibody model using specified parameters and optionally a specified antigenic map and infection history.
#' @param pars the vector of named model parameters, including `boost_long`, `boost_short`,`boost_delay`,`wane_long`,`wane_short`,`cr_long`, and `cr_short`.
#' @param times the vector of times to solve the model over. A continuous vector of discrete timepoints. Can be left to NULL if this information is included in the `antigenic_map` argument.
#' @param infection_history the vector of times matching entries in `times` to simulate infections in.
#' @param antigenic_map the antigenic map to solve the model with. Can be left to NULL to ssume all biomarker IDs have the same antigenic coordinates.
#' @return a data frame with variables `sample_times`, `biomarker_id` and `antibody_level`
#' @examples
#' simulate_antibody_model(c("boost_long"=2,"boost_short"=3,"boost_delay"=1,"wane_short"=0.2,"wane_long"=0.01, "antigenic_seniority"=0,"cr_long"=0.1,"cr_short"=0.03), times=seq(1,25,by=1),infection_history=NULL,antigenic_map=example_antigenic_map)
#'  
#' @export
simulate_antibody_model <- function(pars, 
                                times=NULL, 
                                infection_history=NULL, 
                                antigenic_map=NULL){
  if(is.null(times) & is.null(antigenic_map)){
    stop("Must provide one of times or antigenic_map to give the possible infection times and biomarker IDs over which to solve the model.")
  }
  
  ## If no vector of times provided, take from the antigenic map
  if(is.null(times) & !is.null(antigenic_map)){
    times <- antigenic_map$inf_times
  }
  
  ## If no antigenic map is provided, create a dummy antigenic map where each element has the same antigenic coordinate
  if(is.null(antigenic_map)){
    antigenic_map <- data.frame(x_coord=1,y_coord=1,inf_times=times)
    biomarker_ids <- 0
  } else {
    biomarker_ids <- match(antigenic_map$inf_times[seq_along(times)], antigenic_map$inf_times[seq_along(times)]) -1
  }
  
  ## Setup antigenic map
  use_antigenic_map <- melt_antigenic_coords(antigenic_map[seq_along(times),c("x_coord","y_coord")])
  antigenic_map_long <- matrix(create_cross_reactivity_vector(use_antigenic_map, pars["cr_long"]),ncol=1)
  antigenic_map_short <- matrix(create_cross_reactivity_vector(use_antigenic_map, pars["cr_short"]),ncol=1)
  
  ## If no infection history was provided, setup a dummy infection history vector with only the first entry as an infection
  if(is.null(infection_history)){
    infection_history <- 1
  }
  
  infection_indices <- infection_history-1
  sample_times <- seq(1, max(times),by=1)
  

  start_levels <- rep(0,length(biomarker_ids))
  
  y <- antibody_model_individual_wrapper(pars["boost_long"],pars["boost_short"],pars["boost_delay"],
                                         pars["wane_short"],pars["wane_long"],pars["antigenic_seniority"],
                                         0,
                                         start_levels,
                                         length(times),
                                         infection_history,
                                         infection_indices,
                                         biomarker_ids,
                                         sample_times,
                                         antigenic_map_long,
                                         antigenic_map_short)
  data.frame(sample_times = rep(sample_times,each=length(biomarker_ids)),biomarker_ids = rep(biomarker_ids,length(sample_times)),antibody_level=y)
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
