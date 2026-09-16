# Modified by an AI assistant on 2026-09-16 using GPT-5. Prevented the warning
# for missing expanded starting levels when the user selected the default
# automatic starting-level behaviour.
#'
#' Generate antibody level credible intervals
#'
#' Generates credible intervals on antibody levels and infection histories from an MCMC chain output.
#' @param chain the full MCMC chain to generate antibody level trajectories from, usually `chains$theta_chain`
#' @param infection_histories the MCMC chain for infection histories, usually `chains$inf_chain`
#' @param antibody_data the antibody data frame, with one row per measurement
#' @param demographics optional data frame identifying the demographic group for each individual. This is used when model parameters are stratified by demographics. See the [demographic stratification and covariate vignette](DEMOGRAPHICS_VIGNETTE_LINK).
#' @param individuals the subset of individual IDs to generate credible intervals for
#' @param antigenic_map (optional) a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. The `inf_times` column identifies the circulation or exposure time represented by each map entry. See \code{\link{example_antigenic_map}}
#' @param possible_exposure_times (optional) if no antigenic map is specified, this argument gives the vector of times at which individuals can be infected
#' @param par_tab the model control table specifying the parameters in the MCMC chain
#' @param nsamp number of draws to take from the posterior
#' @param add_residuals if true, returns an extra output summarising residuals between the model prediction and data
#' @param measurement_bias default NULL, optional data frame mapping each `biomarker_id` and `biomarker_group` combination to the `rho_index` of the measurement-shift parameter that it uses. See the [advanced features vignette](ADVANCED_FEATURES_VIGNETTE_LINK).
#' @param for_res_plot TRUE/FALSE value. If using the output of this for plotting of residuals, returns the actual data points rather than summary statistics
#' @param expand_antibody_data TRUE/FALSE value. If TRUE, solves antibody level predictions for every observed biomarker ID at every sample time in the study period. If FALSE, only the biomarker IDs and sample times present in antibody_data are used.
#' @param expand_to_all_times TRUE/FALSE value. If TRUE, uses all possible exposure times as sample times when expanding the prediction data. If FALSE, only the sample times represented in antibody_data are used.
#' @param expand_to_all_biomarker_ids TRUE/FALSE value. If TRUE, solves antibody level predictions for every biomarker ID in the antigenic map while retaining the sample times in antibody_data.
#' @param antibody_level_before_infection TRUE/FALSE value. If TRUE, solves antibody level predictions, but gives the predicted antibody level at a given time point BEFORE any infection during that time occurs.
#' @param for_regression if TRUE, returns posterior draws rather than posterior summaries
#' @param data_type integer, currently accepting 1, 2, or 3. Set to 1 for discrete, bounded data, 2 for continuous, bounded data, or 3 for continuous data with the false-positive observation model. For bounded data, the limits are given by `min_measurement` and `max_measurement` in par_tab.
#' @param start_level `"none"` or a starting-level summary or data frame. A starting level is the antibody level assigned before the modelled infection history begins. With `"none"`, starting levels are set to zero. See the [advanced features vignette](ADVANCED_FEATURES_VIGNETTE_LINK).
#' @param exponential_waning if TRUE, assumes exponential rather than linear waning
#' @return a list with the antibody level predictions (95% credible intervals, median and multivariate posterior mode) and the probabilities of infection for each individual in each epoch
#' @examples
#' \dontrun{
#' data(example_theta_chain)
#' data(example_inf_chain)
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' data(example_par_tab)
#'
#' y <- get_antibody_level_predictions(
#'   chain = example_theta_chain,
#'   infection_histories = example_inf_chain,
#'   antibody_data = example_antibody_data,
#'   individuals = unique(example_antibody_data$individual),
#'   antigenic_map = example_antigenic_map,
#'   par_tab = example_par_tab,
#'   expand_antibody_data = FALSE
#' )
#' }
#' @export
get_antibody_level_predictions <- function(chain, infection_histories, antibody_data,
                                           demographics=NULL,
                                           individuals, antigenic_map=NULL,
                                           possible_exposure_times=NULL, par_tab,
                                           nsamp = 1000, add_residuals = FALSE,
                                           measurement_bias = NULL,
                                           for_res_plot = FALSE, expand_antibody_data = FALSE,
                                           expand_to_all_times=FALSE,
                                           expand_to_all_biomarker_ids=FALSE,
                                           antibody_level_before_infection=FALSE, for_regression=FALSE,
                                           data_type=1,start_level="none",
                                           exponential_waning=FALSE){
  user_supplied_start_levels <- inherits(start_level, c("data.frame", "tibble")) &&
    !isTRUE(attr(start_level, "automatic_start_levels"))
  par_tab <- add_scale_pars(par_tab,antibody_data,demographics)
  ## Get unique demographic groups from full data set, not just the subset
  if(!is.null(demographics)){
    demographic_groups <- create_demographic_table(demographics,par_tab)
  } else {
    demographic_groups <- create_demographic_table(antibody_data,par_tab)
  }
  
  ## Need to align the iterations of the two MCMC chains
  ## and choose some random samples
  samps <- intersect(unique(infection_histories$samp_no), unique(chain$samp_no))
  chain <- chain[chain$samp_no %in% samps, ]
  infection_histories <- infection_histories[infection_histories$samp_no %in% samps, ]
  
  ## Convert samp_no and chainno to a single samp_no index
  if(!("chain_no" %in% colnames(chain))){
    chain$chain_no <- 1
    infection_histories$chain_no <- 1
  }
  chain <- chain %>% dplyr::group_by(chain_no,samp_no) %>% dplyr::mutate(samp_no = cur_group_id()) %>% dplyr::ungroup() %>% dplyr::mutate(chain_no = 1) %>% arrange(samp_no)
  infection_histories <- infection_histories %>% dplyr::group_by(chain_no,samp_no) %>% dplyr::mutate(samp_no = cur_group_id()) %>% dplyr::ungroup() %>% dplyr::mutate(chain_no = 1) %>% arrange(samp_no)
  samps <- unique(chain$samp_no)
  nsamp <- min(nsamp, length(unique(chain$samp_no)))
  
  ## Take subset of individuals
  antibody_data2 <- antibody_data %>% dplyr::select(-c(sample_time,measurement,repeat_number)) %>% distinct()
  antibody_data <- antibody_data[antibody_data$individual %in% individuals, ]
  if(!is.null(demographics)){
    demographics <- demographics[demographics$individual %in% individuals,]
    demographics$individual <- match(demographics$individual, individuals)
  }
  
  infection_histories <- infection_histories[infection_histories$i %in% individuals, ]
  antibody_data$individual <- match(antibody_data$individual, individuals)
  infection_histories$i <- match(infection_histories$i, individuals)
  if(class(start_level) %in% c("data.frame","tibble")){
    start_level$individual <- match(start_level$individual, individuals)
  }
  #if(class(start_level) %in% c("data.frame","tibble")){
  #  start_index_tmp <- start_level$start_index
  #  start_level <- start_level[start_level$individual %in% individuals,]
  #  start_level$individual <- match(start_level$individual, individuals)
  #  start_level$start_index <- 1:nrow(start_level)# match(start_level$start_index, start_index_tmp)
  #}
  ## Format the antigenic map to solve the model 
  ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
  if (!is.null(antigenic_map)) {
    possible_exposure_times_tmp <- unique(antigenic_map$inf_times) 
    ## If possible exposure times was not specified, use antigenic map times instead
    if(is.null(possible_exposure_times)) {
      possible_exposure_times <- possible_exposure_times_tmp
    }
  } else {
    ## Create a dummy map with entries for each observation type
    antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,"inf_times"=possible_exposure_times)
  }
  
  
  nstrain <- length(possible_exposure_times)
  n_indiv <- length(individuals)
  if(!("biomarker_group" %in% colnames(antibody_data))){
    antibody_data$biomarker_group <- 1
  }
  unique_biomarker_groups <- unique(antibody_data$biomarker_group)
  
  ## Empty data structures to save output to
  infection_history_dens <- NULL
  tmp_samp <- sample(samps, nsamp)
  
  ## See the function in posteriors.R
  antibody_data1 <- antibody_data
  start_level1 <- start_level
  if (expand_antibody_data | expand_to_all_biomarker_ids) {
    if(expand_to_all_times){
      expand_times <- possible_exposure_times
    } else {
      expand_times <- unique(antibody_data$sample_time)
    }
    if(expand_to_all_biomarker_ids){
      if("biomarker_group" %in% colnames(antigenic_map)){
        biomarker_data <- antigenic_map %>%
          dplyr::transmute(biomarker_group, biomarker_id = inf_times) %>%
          dplyr::distinct()
      } else {
        biomarker_data <- tidyr::expand_grid(
          biomarker_group = unique(antibody_data$biomarker_group),
          biomarker_id = unique(antigenic_map$inf_times)
        )
      }
      antibody_data1 <- expand.grid(
        individual = unique(antibody_data$individual),
        sample_time = expand_times,
        biomarker_group=unique(antibody_data$biomarker_group),
        measurement = 0, repeat_number = 1
      ) %>%
        dplyr::left_join(biomarker_data, by="biomarker_group") %>%
        dplyr::left_join(
          antibody_data %>%
            dplyr::select(-c(measurement,repeat_number,biomarker_id,biomarker_group)) %>%
            dplyr::distinct(),
          by=c("individual","sample_time")
        )
    } else {
      antibody_data1 <- expand.grid(
        individual = unique(antibody_data$individual),
        sample_time = expand_times,
        biomarker_group=unique(antibody_data$biomarker_group),
        measurement = 0, repeat_number = 1
      )
      antibody_data2 <- antibody_data %>% dplyr::select(-c(sample_time,measurement,repeat_number)) %>% distinct()
      antibody_data1 <- merge(antibody_data1, antibody_data2)
    }
    antibody_data1 <- antibody_data1 %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
    
    ## Create full start level data
    if(expand_to_all_biomarker_ids){
      if(class(start_level) %in% c("data.frame","tibble")){
        start_level_complete <- start_level %>%
          dplyr::select(individual, biomarker_id, biomarker_group, starting_level, start_index) %>%
          dplyr::distinct()
      } else {
        start_level_summary <- if(is.character(start_level)) start_level else "none"
        start_level_complete <- create_start_level_data(antibody_data, start_level_summary, FALSE) %>%
          dplyr::select(individual, biomarker_id, biomarker_group, starting_level, start_index) %>%
          dplyr::distinct()
      }
      missing_start_levels <- antibody_data1 %>%
        dplyr::select(individual, biomarker_id, biomarker_group) %>%
        dplyr::distinct() %>%
        dplyr::anti_join(start_level_complete,
                         by=c("individual","biomarker_id","biomarker_group"))
      # Warn only when a user-supplied table is incomplete; automatic summaries use zero silently.
      if(nrow(missing_start_levels) > 0 && user_supplied_start_levels){
        warning(paste0(
          "No starting levels were supplied for ", nrow(missing_start_levels),
          " individual-biomarker combinations created by expanding to all biomarker IDs. ",
          "These starting levels have been set to zero. This may give unexpected predictions,",
          " particularly when an unobserved biomarker ID lies between observed IDs with non-zero starting levels."
        ), call.=FALSE)
      }
      if(nrow(missing_start_levels) > 0){
        missing_start_levels <- missing_start_levels %>%
          dplyr::mutate(starting_level=0,
                        start_index=max(start_level_complete$start_index) + dplyr::row_number())
        start_level_complete <- dplyr::bind_rows(start_level_complete, missing_start_levels)
      }
      start_level1 <- antibody_data1 %>%
        dplyr::left_join(start_level_complete,
                         by=c("individual","biomarker_id","biomarker_group"))
    } else {
      start_level1 <- antibody_data1 %>% left_join(start_level %>% select(individual,biomarker_id,biomarker_group,starting_level, start_index) %>% distinct(),
                                                   by=c("individual","biomarker_group","biomarker_id"))
    }
    
  }
  antibody_data1 <- antibody_data1 %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
  
  model_func <- create_posterior_func(par_tab, antibody_data1, antigenic_map, possible_exposure_times,
                                      prior_version=2,
                                      measurement_bias = measurement_bias, function_type = 4,
                                      antibody_level_before_infection=antibody_level_before_infection,
                                      data_type=data_type,start_level=start_level1,
                                      demographics=demographics,demographic_groups=demographic_groups,
                                      exponential_waning=exponential_waning
  )
  
  predicted_titres <- residuals <- residuals_floor <- 
    observed_predicted_titres <- matrix(nrow = nrow(antibody_data1), ncol = nsamp)
  samp_record <- numeric(nsamp)
  
  
  ## For each sample, take values for theta and infection histories and simulate titres
  inf_hist_all <- list(nsamp)
  for (i in 1:nsamp) {
    index <- tmp_samp[i]
    pars <- get_index_pars(chain, samp_no=index,chain_no=1)
    pars <- pars[!(names(pars) %in% c("posterior_prob", "likelihood", "prior_prob",
                                      "samp_no", "total_infections", "chain_no"
    ))]
    names(pars) <- par_tab$names
    tmp_inf_hist <- infection_histories[infection_histories$samp_no == index, ]
    tmp_inf_hist <- as.matrix(Matrix::sparseMatrix(i = tmp_inf_hist$i, j = tmp_inf_hist$j, x = tmp_inf_hist$x, dims = c(n_indiv, nstrain)))
    predicted_titres[, i] <- model_func(pars, tmp_inf_hist)
    for(biomarker_group in unique_biomarker_groups){
      observed_predicted_titres[which(antibody_data1$biomarker_group == biomarker_group),i] <- add_noise(predicted_titres[which(antibody_data1$biomarker_group == biomarker_group),i], pars, NULL, NULL,data_type=data_type[biomarker_group])
    }
    inf_hist_all[[i]] <- tmp_inf_hist
    ## Get residuals between observations and predictions
    residuals[, i] <- antibody_data1$measurement - predicted_titres[, i]
    residuals_floor[,i] <- antibody_data1$measurement - observed_predicted_titres[,i]
    samp_record[i] <- index
  }
  
  colnames(predicted_titres) <- tmp_samp
  ## If generating for residual plot, can return now
  if (for_res_plot) {
    return(list(residuals, samp_record, antibody_data1,
                predicted_titres,
                observed_predicted_titres,
                residuals_floor))
  }
  
  ## Get 95% credible interval and means
  dat2 <- t(apply(predicted_titres, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))))
  
  ## Get 95% credible interval and means of observations
  obs_dat <- t(apply(observed_predicted_titres, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))))
  
  residuals <- t(apply(residuals, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.5, 0.975)))))
  residuals <- cbind(antibody_data1, residuals)
  
  residuals_obs <- t(apply(residuals_floor, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.5, 0.975)))))
  residuals_obs <- cbind(antibody_data1, residuals_obs)
  
  ## Find multivariate posterior mode estimate from the chain
  best_pars <- get_best_pars(chain)
  best_pars <- best_pars[!(names(best_pars) %in% c(
    "posterior_prob", "likelihood", "prior_prob",
    "samp_no", "total_infections", "chain_no"
  ))]
  #best_pars <- best_pars[names(best_pars) %in% par_tab$names]
  best_I <- chain$samp_no[which.max(chain$posterior_prob)]
  best_inf <- infection_histories[infection_histories$samp_no == best_I, ]
  best_inf <- as.matrix(Matrix::sparseMatrix(i = best_inf$i, j = best_inf$j, x = best_inf$x, dims = c(n_indiv, nstrain)))
  ## Generate trajectory for best parameters
  best_traj <- model_func(best_pars, best_inf)
  best_residuals <- antibody_data1$measurement - best_traj
  best_residuals <- cbind(antibody_data1, best_residuals, "samp_no" = best_I)
  dat2 <- as.data.frame(dat2)
  obs_dat <- as.data.frame(obs_dat)
  
  colnames(dat2) <- colnames(obs_dat) <- c("mean","lower", "lower_50", "median", "upper_50", "upper")
  dat2$max <- best_traj
  dat2 <- cbind(antibody_data1, dat2)
  obs_dat <- cbind(antibody_data1, obs_dat)
  tmp_inf_chain <- data.table(subset(infection_histories, samp_no %in% tmp_samp))
  
  ## Get infection history density for each individual and each epoch
  data.table::setkey(tmp_inf_chain, "i", "j")
  infection_history_dens <- tmp_inf_chain[, list(V1 = sum(x) / length(tmp_samp)), by = key(tmp_inf_chain)]
  infection_history_dens$j <- possible_exposure_times[infection_history_dens$j]
  colnames(infection_history_dens) <- c("individual", "variable", "value")
  infection_history_final <- infection_history_dens
  best_inf <- data.frame(best_inf)
  best_inf$individual <- 1:nrow(best_inf)
  best_inf$individual <- individuals[best_inf$individual]
  
  dat2$individual <- individuals[dat2$individual]
  infection_history_final$individual <- individuals[infection_history_final$individual]
  if(for_regression){
    return(list("all_predictions"=predicted_titres, "all_predictions_obs"=observed_predicted_titres,
                "all_inf_hist"=inf_hist_all,
                "summary_titres"=dat2,"best_inf_hist"=best_inf, "predicted_observations"=obs_dat)) 
  }
  
  if (add_residuals) {
    result <- list("predictions" = dat2, "histories" = infection_history_final, 
                   "residuals" = residuals, "bestRes" = best_residuals,"best_infhist"=best_inf,
                   "predicted_observations"=obs_dat)
  } else {
    result <- list("predictions" = dat2, "histories" = infection_history_final,
                   "best_infhist"=best_inf, "predicted_observations"=obs_dat)
  }
  return(result)
}
