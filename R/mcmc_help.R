#' Generate starting parameter table
#'
#' Generates a version of \code{par_tab} with random values between \code{lower_start} and \code{upper_start}
#' @param par_tab See \code{\link{example_par_tab}}
#' @return a data frame matching par_tab
#' @family mcmc
#' @examples
#' data(example_par_tab)
#' start_tab <- generate_start_tab(example_par_tab)
#' @export
generate_start_tab <- function(par_tab){
    for(i in 1:nrow(par_tab)){
        if(par_tab[i,"fixed"] == 0){
            par_tab[i, "values"] <- runif(1,par_tab[i,"lower_start"], par_tab[i, "upper_start"])
        }
    }
    return(par_tab)        
}


#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param target_acceptance_rate_theta the desired acceptance rate
#' @param pcur the current acceptance rate
#' @return the scaled step size
#' @export
#' @family mcmc
#' @useDynLib serosolver
scaletuning <- function(step, target_acceptance_rate_theta, pcur) {
  if(is.finite(pcur)){
    if (pcur == 1) pcur <- 0.99
    if (pcur == 0) pcur <- 0.01
    step <- (step * qnorm(target_acceptance_rate_theta / 2)) / qnorm(pcur / 2)
    if (step > 1) step <- 1
    step <- max(0.00001, step)
  }
  return(step)
}

#' Robins and Monro scaler, thanks to Michael White
#' @family mcmc
#' @export
rm_scale <- function(step_scale, mc, target_acceptance_rate_theta, log_prob, N_adapt) {
  dd <- exp(log_prob)
  if (dd < -30) {
    dd <- 0
  }
  dd <- min(dd, 1)

  rm_temp <- (dd - target_acceptance_rate_theta) / ((mc + 1) / (0.01 * N_adapt + 1))

  out <- step_scale * exp(rm_temp)

  out <- max(out, 0.02)
  out <- min(out, 2)
  out
}

#' Generate random infection history matrix
#'
#' Creates a random infection history matrix by sampling from the infection history model under prior version 2
#' @param antibody_data the data frame of antibody data. See \code{\link{example_antibody_data}}
#' @param possible_exposure_times the vector of times that individuals can be infected
#' @param infection_model_prior_shape1 shape parameter 1 (alpha) of the Beta distribution
#' @param infection_model_prior_shape2 shape parameter 2 (beta) of the Beta distribution
#' @return an n (number of individuals) by m (number of possible exposure times) matrix containing 1s and 0s, representing infections
#' @family setup_infection_histories
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_times
#' setup_infection_histories_prior(example_antibody_data, times, 1, 1)
#' @export
setup_infection_histories_prior <- function(antibody_data, possible_exposure_times, infection_model_prior_shape1=1,infection_model_prior_shape2=1){
  DOBs <- unique(antibody_data[, c("individual", "birth")])[, 2]
  n_indiv <- length(unique(antibody_data$individual))
  n_infections <- length(possible_exposure_times)
  
  age_mask <- create_age_mask(DOBs, possible_exposure_times)
  sample_mask <- create_sample_mask(antibody_data, possible_exposure_times)
  masks <- data.frame(cbind(age_mask, sample_mask))
  
  ## Go through each possible infection and assign random with p total_infs/n_alive
  infection_histories <- matrix(0, nrow = n_indiv, ncol = n_infections)
  
  probs <- rbeta(length(possible_exposure_times),infection_model_prior_shape1,infection_model_prior_shape2)
  
  for (i in 1:nrow(masks)) {
    years <- age_mask[i]:sample_mask[i]
    n <- length(years)
    infection_histories[i, years] <- rbinom(length(years), 1,probs[years])
  }
  return(infection_histories)
}

#' Propose initial infection histories based on antibody levels
#'
#' Very similar to \code{\link{setup_infection_histories}}, but is not restricted to placing starting infections against antigens to which an individual has a measurable antibody level Given a matrix of antibody data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' The idea is to move along time and look at an individual's antibody level against each antigen/variant Where antibodies are elevated, this suggests an infection. However, to avoid suggesting multiple infections for regions of high antigenic similarity, we place a necessary gap (defined by `space`) between proposed infection times.
#' @param antibody_data the matrix of titres data with columns for individual, sample, and antibody level
#' @param possible_exposure_times vector of real times for all strains
#' @param space how many epochs must separate proposed infections
#' @param antibody_cutoff specifies how high the antibody level must be to imply an infection
#' @param sample_prob if antibody levels suggest an infection, then add an infection with 1 minus this probability
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of time periods for potential infection
#' @family setup_infection_histories
#' @examples
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' start_inf <- setup_infection_histories_antibody_level(example_antibody_data, example_antigenic_map$inf_times)
#' @export
setup_infection_histories_antibody_level <- function(antibody_data, possible_exposure_times, space = 5, antibody_cutoff = 2, sample_prob = 0.9) {
  start_inf <- NULL
  individuals <- unique(antibody_data$individual)
  ages <- unique(antibody_data[, c("individual", "birth")])
  
  ## For each individual
  for (individual in individuals) {
    ## Isolate that individual's data and date of birth
    dat <- antibody_data[antibody_data$individual == individual, ]
    unique_biomarkers <- unique(dat$biomarker_id)
    dob <- as.numeric(ages[ages$individual == individual, "birth"])
    times <- possible_exposure_times

    ## What was the most recent strain that the individual could get
    sample_mask <- create_sample_mask(dat, possible_exposure_times)
    
    ## Only look at biomarker IDs that circulated when an individual was alive and for samples not in the future
    times <- times[times >= dob & times <= possible_exposure_times[sample_mask]]
    
    inf_times <- NULL
    i <- 0
    while (i < length(times)) { ## Start going through each biomarker_id
      i <- i + 1
      biomarker_id <- times[i] ## Get current biomarker_id of interest
      dist <- 0
      measured_biomarker <- unique_biomarkers[which(abs(unique_biomarkers - biomarker_id) == min(abs(unique_biomarkers - biomarker_id)))][1]
      measurement <- max(dat[dat$biomarker_id == measured_biomarker, "measurement"]) ## Get max antibody level against this biomarker_id
      if (measurement >= antibody_cutoff) { ## If elevated against this biomarker_id, assume an infection
        new_inf <- biomarker_id
        ## Begin counting up distance
        while (dist < space & i < length(times)) {
          i <- i + 1
          biomarker_id <- times[i]
          dist <- biomarker_id - new_inf ## How many years since last infection?
          measured_biomarker <- unique_biomarkers[which(abs(unique_biomarkers - biomarker_id) == min(abs(unique_biomarkers - biomarker_id)))][1]
          new_measurement <- max(dat[dat$biomarker_id == measured_biomarker, "measurement"]) ## Get max measurement against next biomarker ID along
          ## If this one is better, replace and reset distance
          if (new_measurement > measurement) {
            new_inf <- measured_biomarker
            measurement <- new_measurement
            dist <- 0
          }
        }
        if (runif(1) > sample_prob) inf_times <- c(inf_times, new_inf)
        dist <- 0
      }
    }
    infections <- rep(0, length(possible_exposure_times))
    infections[match(inf_times, possible_exposure_times)] <- 1
    start_inf <- rbind(start_inf, infections)
  }
  colnames(start_inf) <- possible_exposure_times
  rownames(start_inf) <- NULL
  return(start_inf)
}


#' Propose initial infection histories
#'
#' Given a matrix of antibody data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' The idea is to move along time in the context of antigenic drift and look at an individual's antibody level against each biomarker ID Where antibody levels are raised, we suggest an infection. However, to avoid suggesting multiple infections for regions of high antigenic similarity, we place a necessary gap (defined by `space`) between proposed infection times.
#' @inheritParams setup_infection_histories_antibody_level
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of time periods for potential infection
#' @family setup_infection_histories
#' @examples
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' start_inf <- setup_infection_histories(example_antibody_data, example_antigenic_map$inf_times)
#' @export
setup_infection_histories <- function(antibody_data, possible_exposure_times, space = 5, antibody_cutoff = 2, sample_prob = 0.9) {
  start_inf <- NULL
  individuals <- unique(antibody_data$individual)
  ages <- unique(antibody_data[, c("individual", "birth")])

  ## For each individual
  for (individual in individuals) {
    ## Isolate that individual's data and date of birth
    dat <- antibody_data[antibody_data$individual == individual, ]
    dob <- as.numeric(ages[ages$individual == individual, "birth"])
    biomarker_ids <- unique(dat$biomarker_id)

    ## What was the most recent infection that the individual could get
    sample_mask <- create_sample_mask(dat, possible_exposure_times)

    ## Only look at biomarker IDs that circulated when an individual was alive and for samples not in the future
    biomarker_ids <- biomarker_ids[biomarker_ids >= dob & biomarker_ids <= possible_exposure_times[sample_mask]]

    inf_times <- NULL
    i <- 0
    while (i < length(biomarker_ids)) { ## Start going through each biomarker_id
      i <- i + 1
      biomarker_id <- biomarker_ids[i] ## Get current biomarker_id of interest
      dist <- 0
      antibody_level <- max(dat[dat$biomarker_id == biomarker_id, "measurement"]) ## Get max measurement against this biomarker_id
      if (antibody_level >= antibody_cutoff) { ## If elevated against this biomarker_id, assume an infection
        new_inf <- biomarker_id
        ## Begin counting up distance
        while (dist < space & i < length(biomarker_ids)) {
          i <- i + 1
          biomarker_id <- biomarker_ids[i]
          dist <- biomarker_id - new_inf ## How many years since last infection?
          new_measurement <- max(dat[dat$biomarker_id == biomarker_id, "measurement"]) ## Get max measurement against next biomarker_id along
          ## If this one is better, replace and reset distance
          if (new_measurement > antibody_level) {
            new_measurement <- biomarker_id
            titre <- new_measurement
            dist <- 0
          }
        }
        if (runif(1) > sample_prob) inf_times <- c(inf_times, new_inf)
        dist <- 0
      }
    }
    infections <- rep(0, length(possible_exposure_times))
    infections[match(inf_times, possible_exposure_times)] <- 1
    start_inf <- rbind(start_inf, infections)
  }
  colnames(start_inf) <- possible_exposure_times
  rownames(start_inf) <- NULL
  return(start_inf)
}




#' Write given infection history to disk
#'
#' @param infection_history the infection history matrix
#' @param file the file location to save to
#' @param samp_no which sample number is this matrix?
#' @param append if TRUE, just adds to the bottom of the file
#' @param col_names if TRUE, saves column names first (only set to true if append = FALSE)
#' @return nothing
#' @family mcmc
#' @export
save_infection_history_to_disk <- function(infection_history, file, samp_no, append = TRUE, col_names = FALSE) {
  save_inf_hist <- Matrix::Matrix(infection_history, sparse = TRUE)
  save_inf_hist <- as.data.frame(Matrix::summary(save_inf_hist))
  if (nrow(save_inf_hist) > 0) {
    save_inf_hist$samp_no <- samp_no
    try(data.table::fwrite(save_inf_hist, file = file, col.names = col_names, row.names = FALSE, sep = ",", append = append))
  }
}

#' Expand sparse infection history matrix
#'
#' @param inf_chain the data table with the saved sparse infection history matrix
#' @param j optional vector of js to expand the infection history chain for
#' @return long format, full infection history matrix chain
#' @export
expand_summary_inf_chain <- function(inf_chain, j_vec = NULL) {
  if (is.null(j_vec)) j_vec <- 1:max(inf_chain$j)
  full_inf_chain <- data.table::CJ(i = min(inf_chain$i):max(inf_chain$i), j = j_vec, samp_no = sort(unique(inf_chain$samp_no)))
  inf_chain <- data.table::data.table(apply(inf_chain, 2, as.numeric))
  summary_with_non_infections <- merge(inf_chain, full_inf_chain, by = c("samp_no", "j", "i"), all = TRUE)
  summary_with_non_infections[is.na(summary_with_non_infections$x), "x"] <- 0
  colnames(summary_with_non_infections) <- c("samp_no", "j", "individual", "x")
  expanded_chain <- data.table::dcast(summary_with_non_infections, samp_no + individual ~ j, value.var = "x")
  return(expanded_chain)
}
