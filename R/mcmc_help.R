#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param popt the desired acceptance rate
#' @param pcur the current acceptance rate
#' @return the scaled step size
#' @export
#' @useDynLib serosolver
scaletuning <- function(step, popt, pcur) {
  if (pcur == 1) pcur <- 0.99
  if (pcur == 0) pcur <- 0.01
  step <- (step * qnorm(popt / 2)) / qnorm(pcur / 2)
  if (step > 1) step <- 1
  step <- max(0.00001, step)
  return(step)
}

#' Robins and Monro scaler, thanks to Michael White
#' @export
rm_scale <- function(step_scale, mc, popt, log_prob, N_adapt) {
  dd <- exp(log_prob)
  if (dd < -30) {
    dd <- 0
  }
  dd <- min(dd, 1)

  rm_temp <- (dd - popt) / ((mc + 1) / (0.01 * N_adapt + 1))

  out <- step_scale * exp(rm_temp)

  out <- max(out, 0.02)
  out <- min(out, 2)
  out
}


#' Propose initial infection histories - OLD VERSION
#'
#' Given a matrix of titre data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' NOTE - MIGHT NEED TO UPDATE THIS FOR GROUPS
#' @param titre_dat the matrix of titres data with columns for individual, sample, and titre
#' @param strain_isolation_times vector of real times for all strains
#' @param sample_prob given an infection seems likely based on titre, suggest infection with 1 minus this probability
#' @param titre_cutoff specifies how high the titre must be to imply an infection
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of potential infecting strains
#' @family setup_infection_histories
#' @export
setup_infection_histories_old <- function(titre_dat, strain_isolation_times, sample_prob, titre_cutoff = 3) {
  SAMPLE_PROB <- sample_prob
  n_indiv <- length(unique(titre_dat$individual))
  n_strain <- length(strain_isolation_times)
  sampling_times <- unique(titre_dat$samples)
  sample_time <- max(sampling_times)
  infection_histories <- matrix(0, nrow = n_indiv, ncol = n_strain)
  strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)

  ages <- unique(titre_dat[, c("individual", "DOB")])
  age_mask <- create_age_mask(ages$DOB, strain_isolation_times)

  strains <- unique(titre_dat$virus)

  index <- 1
  ## For each individual
  for (indiv in unique(titre_dat$individual)) {
    ## For each sampling time
    tmp_inf_hist <- numeric(n_strain)
    index2 <- 1
    tmp_strains <- strains[age_mask[indiv]:strain_mask[indiv]]

    for (strain in tmp_strains) {
      tmp_titre <- max(titre_dat[titre_dat$virus == strain &
        titre_dat$individual == indiv, "titre"], 0)
      tmp_inf <- 0
      # If high titre, set infection presence to 1 with some probability (0.2)
      if (!is.na(tmp_titre) && (tmp_titre >= titre_cutoff & runif(1) > SAMPLE_PROB)) {
        tmp_inf <- 1
      }
      tmp_inf_hist[index2] <- tmp_inf
      index2 <- index2 + 1
    }
    infection_histories[index, ] <- tmp_inf_hist
    ## Add infection at some point in the last 10 years
    forced_infection <- which(tmp_strains == sample(tmp_strains[tmp_strains <= sample_time & tmp_strains >= (sample_time - 10)], 1))
    ## Pick strain within plausible region to add
    if (sum(infection_histories[index, which(tmp_strains < sample_time)]) == 0) infection_histories[index, forced_infection] <- 1

    index <- index + 1
  }
  infection_histories
}



#' Initial infection history prior on total
#'
#' Sets up an initial random infection history, where start is sampled from a prior on the total number of infections across all individuals and years
#' @inheritParams setup_infection_histories_new_2
#' @param alpha alpha parameter for beta distribution to sample from
#' @param beta beta parameter for beta distribution to sample from
#' @return an infection history matrix
#' @family setup_infection_histories
#' @export
setup_infection_histories_total <- function(titre_dat, strain_isolation_times, alpha = 1, beta = 1) {
  DOBs <- unique(titre_dat[, c("individual", "DOB")])[, 2]
  n_indiv <- length(unique(titre_dat$individual))
  n_strain <- length(strain_isolation_times)

  age_mask <- create_age_mask(DOBs, strain_isolation_times)
  strain_mask <- create_strain_mask(titre_dat, strain_isolation_times)
  masks <- data.frame(cbind(age_mask, strain_mask))
  ## Number of people that were born before each year and have had a sample taken since that year happened
  n_alive <- sum(sapply(seq(1, length(strain_isolation_times)), function(x) nrow(masks[masks$age_mask <= x & masks$strain_mask >= x, ])))

  ## How many infections are we spreading across everyone?
  total_infs <- rbb(1, n_alive, alpha, beta)

  ## Go through each possible infection and assign random with p total_infs/n_alive
  infection_histories <- matrix(0, nrow = n_indiv, ncol = n_strain)

  for (i in 1:nrow(masks)) {
    years <- age_mask[i]:strain_mask[i]
    n <- length(years)
    infection_histories[i, years] <- ifelse(runif(n) < total_infs / n_alive, 1, 0)
  }
  return(infection_histories)
}


#' Propose initial infection histories
#'
#' Given a matrix of titre data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' The idea is to move along time in the context of antigenic drift and look at an individual's titre against each strain. Where titres are raised, we suggest an infection. However, to avoid suggesting multiple infections for regions of high antigenic similarity, we place a necessary gap (defined by `space`) between proposed infection times.
#' @inheritParams setup_infection_histories_new_2 
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of potential infecting strains
#' @export
setup_infection_histories_new <- function(titre_dat, strain_isolation_times, space = 5, titre_cutoff = 2, sample_prob=0.9) {
  start_inf <- NULL
  individuals <- unique(titre_dat$individual)
  ages <- unique(titre_dat[, c("individual", "DOB")])

  ## For each individual
  for (individual in individuals) {
    ## Isolate that individual's data and date of birth
    dat <- titre_dat[titre_dat$individual == individual, ]
    dob <- as.numeric(ages[ages$individual == individual, "DOB"])
    strains <- unique(dat$virus)
    # strains <- strain_isolation_times
    ## What was the most recent strain that the individual could get
    strain_mask <- create_strain_mask(dat, strain_isolation_times)

    ## Only look at strains that circulated when an individual was alive and for samples not in the future
    strains <- strains[strains >= dob & strains <= strain_isolation_times[strain_mask]]

    inf_years <- NULL
    i <- 0
    while (i < length(strains)) { ## Start going through each strain
      i <- i + 1
      strain <- strains[i] ## Get current strain of interest
      dist <- 0
      titre <- max(dat[dat$virus == strain, "titre"]) ## Get max titre against this strain
      if (titre >= titre_cutoff) { ## If elevated against this strain, assume an infection
        new_inf <- strain
        ## Begin counting up distance
        while (dist < space & i < length(strains)) {
          i <- i + 1
          strain <- strains[i]
          dist <- strain - new_inf ## How many years since last infection?
          new_titre <- max(dat[dat$virus == strain, "titre"]) ## Get max titre against next strain along
          ## If this one is better, replace and reset distance
          if (new_titre > titre) {
            new_inf <- strain
            titre <- new_titre
            dist <- 0
          }
        }
        if(runif(1) > sample_prob) inf_years <- c(inf_years, new_inf)
        dist <- 0
      }
    }
    infections <- rep(0, length(strain_isolation_times))
    infections[match(inf_years, strain_isolation_times)] <- 1
    start_inf <- rbind(start_inf, infections)
  }
  colnames(start_inf) <- strain_isolation_times
  rownames(start_inf) <- NULL
  return(start_inf)
}

#' Propose initial infection histories 2 - use this!
#'
#' Very similar to \code{\link{setup_infection_histories_new}}, but is not restricted to placing starting infections against viruses to which an individual has a titre. Given a matrix of titre data, proposes plausible initial infection histories from which to begin MCMC sampling.
#' The idea is to move along time in the context of antigenic drift and look at an individual's titre against each strain. Where titres are raised, we suggest an infection. However, to avoid suggesting multiple infections for regions of high antigenic similarity, we place a necessary gap (defined by `space`) between proposed infection times.
#' @param titre_dat the matrix of titres data with columns for individual, sample, and titre
#' @param strain_isolation_times vector of real times for all strains
#' @param space how many epochs must separate proposed infections
#' @param titre_cutoff specifies how high the titre must be to imply an infection
#' @param sample_prob if titre suggests an infection, then add an infection with 1 minus this probability
#' @return an nxm matrix of infection histories containing 1s and 0s, where n is the number of individuals and m is the number of potential infecting strains
#' @family setup_infection_histories
#' @export
setup_infection_histories_new_2 <- function(titre_dat, strain_isolation_times, space = 5, titre_cutoff = 2, sample_prob=0.9) {
  start_inf <- NULL
  individuals <- unique(titre_dat$individual)
  ages <- unique(titre_dat[, c("individual", "DOB")])

  ## For each individual
  for (individual in individuals) {
    ## Isolate that individual's data and date of birth
    dat <- titre_dat[titre_dat$individual == individual, ]
    unique_viruses <- unique(dat$virus)
    dob <- as.numeric(ages[ages$individual == individual, "DOB"])
    strains <- strain_isolation_times
    # strains <- strain_isolation_times
    ## What was the most recent strain that the individual could get
    strain_mask <- create_strain_mask(dat, strain_isolation_times)

    ## Only look at strains that circulated when an individual was alive and for samples not in the future
    strains <- strains[strains >= dob & strains <= strain_isolation_times[strain_mask]]

    inf_years <- NULL
    i <- 0
    while (i < length(strains)) { ## Start going through each strain
      i <- i + 1
      strain <- strains[i] ## Get current strain of interest
      dist <- 0
      measured_virus <- unique_viruses[which(abs(unique_viruses - strain) == min(abs(unique_viruses - strain)))][1]
      titre <- max(dat[dat$virus == measured_virus, "titre"]) ## Get max titre against this strain
      if (titre >= titre_cutoff) { ## If elevated against this strain, assume an infection
        new_inf <- strain
        ## Begin counting up distance
        while (dist < space & i < length(strains)) {
          i <- i + 1
          strain <- strains[i]
          dist <- strain - new_inf ## How many years since last infection?
          measured_virus <- unique_viruses[which(abs(unique_viruses - strain) == min(abs(unique_viruses - strain)))][1]
          new_titre <- max(dat[dat$virus == measured_virus, "titre"]) ## Get max titre against next strain along
          ## If this one is better, replace and reset distance
          if (new_titre > titre) {
            new_inf <- strain
            titre <- new_titre
            dist <- 0
          }
        }
        if(runif(1) > sample_prob) inf_years <- c(inf_years, new_inf)
        dist <- 0
      }
    }
    infections <- rep(0, length(strain_isolation_times))
    infections[match(inf_years, strain_isolation_times)] <- 1
    start_inf <- rbind(start_inf, infections)
  }
  colnames(start_inf) <- strain_isolation_times
  rownames(start_inf) <- NULL
  return(start_inf)
}


#' Write given infection history to disk
#'
#' @param infection_history the infection history matrix
#' @param file the file location to save to
#' @param sampno which sample number is this matrix?
#' @param append if TRUE, just adds to the bottom of the file
#' @param col_names if TRUE, saves column names first (only set to true if append = FALSE)
#' @return nothing
#' @export
save_infection_history_to_disk <- function(infection_history, file, sampno, append = TRUE, col_names = FALSE) {
  save_inf_hist <- Matrix::Matrix(infection_history, sparse = TRUE)
  save_inf_hist <- as.data.frame(Matrix::summary(save_inf_hist))
  if (nrow(save_inf_hist) > 0) {
      save_inf_hist$sampno <- sampno
      try(data.table::fwrite(save_inf_hist, file = file, col.names = col_names, row.names = FALSE, sep = ",", append = append))
  }
}

#' Expand sparse infection history matrix
#'
#' @param inf_chain the data table with the saved sparse infection history matrix
#' @param j optional vector of js to expand the infection history chain for
#' @return long format, full infection history matrix chain
#' @export
expand_summary_infChain <- function(inf_chain,j_vec=NULL){
    if(is.null(j_vec)) j_vec<-1:max(inf_chain$j)
    full_inf_chain <- data.table::CJ(i=min(inf_chain$i):max(inf_chain$i), j=j_vec, sampno=sort(unique(inf_chain$sampno)))
    inf_chain <- data.table::data.table(apply(inf_chain, 2, as.numeric))
    summary_with_non_infections <- merge(inf_chain,full_inf_chain,by=c("sampno","j","i"),all=TRUE)
    summary_with_non_infections[is.na(summary_with_non_infections$x),"x"] <- 0
    colnames(summary_with_non_infections) <- c("sampno","j","individual","x")
    expanded_chain <- data.table::dcast(summary_with_non_infections, sampno + individual ~ j,value.var="x")
    return(expanded_chain)
}
