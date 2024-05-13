#' Read in start_tab
#'
#' Searches the working directory for a file called "start_tab.csv" and reads in the first one found
#' @param location string giving the relative file path to be searched
#' @return data frame, the first found starting parameter tables in the current working directory
#' @family load_data_functions
#' @examples
#' \dontrun{load_start_tab()}
#' @export
load_start_tab <- function(location = getwd()) {
  files <- Sys.glob(file.path(location, "*_start_tab.csv"))
  files_old <- Sys.glob(file.path(location, "*_startTab.csv"))
  files <- c(files, files_old)
  if (length(files) < 1) {
    message("Error - no files found")
    return(NULL)
  }
  start_tab <- read.csv(files[1], stringsAsFactors = FALSE)
  return(start_tab)
}

#' Read in antibody_data
#'
#' Searches the working directory for a file with "antibody_data.csv" in the name and reads in the first one found
#' @inheritParams load_start_tab
#' @return data frame, the first found starting parameter tables in the current working directory
#' @family load_data_functions
#' @examples
#' \dontrun{load_antibody_data()}
#' @export
load_antibody_data <- function(location = getwd()) {
  files <- Sys.glob(file.path(location, "*_antibody_data.csv"))
  files_antibody_data <- Sys.glob(file.path(location, "*_titre_dat.csv"))
  files_titreDat <- Sys.glob(file.path(location, "*_titreDat.csv"))
  files <- c(files, files_antibody_data, files_titreDat)
  if (length(files) < 1) {
    message("Error - no files found")
    return(NULL)
  }
  antibody_data <- read.csv(files[1], stringsAsFactors = FALSE)
  return(antibody_data)
}


#' Read in antigenic_map
#'
#' Searches the working directory for a file with "antigenic_map.csv" in the name and reads in the first one found
#' @inheritParams load_start_tab
#' @return data frame, the first found antigenic map tables in the current working directory
#' @family load_data_functions
#' @examples
#' \dontrun{load_antigenic_map_file()}
#' @export
load_antigenic_map_file <- function(location = getwd()) {
  files <- Sys.glob(file.path(location, "*_antigenic_map.csv"))
  files_old <- Sys.glob(file.path(location, "*_antigenicMap.csv"))
  files <- c(files, files_old)
  if (length(files) < 1) {
    message("Error - no files found")
    return(NULL)
  }
  antigenic_map <- read.csv(files[1], stringsAsFactors = FALSE)
  return(antigenic_map)
}




#' Combine theta and infection history chains
#'
#' Reads in all MCMC chains for theta and infection histories from the specified directory, adding in the total number of infections
#' @inheritParams load_theta_chains
#' @return a list of the concatenated and individual chains (4 elements, either data frames of coda::mcmc objects)
#' @family load_data_functions
#' @examples
#' \dontrun{load_mcmc_chains(par_tab=par_tab, estimated_only=TRUE,thin=10,burnin=5000,convert_mcmc=TRUE)}
#' @export
load_mcmc_chains <- function(location = getwd(), par_tab = NULL, estimated_only = FALSE, thin = 1, burnin = 0, convert_mcmc = FALSE) {
    ## Load in theta chains
    theta_chains <- load_theta_chains(location, par_tab, estimated_only, thin, burnin)
    ## Load in infection history chains
    inf_chains <- load_infection_chains(location, thin, burnin)

    chain <- theta_chains$chain
    inf_chain <- inf_chains$chain
    theta_list_chains <- theta_chains$list
    inf_list_chains <- inf_chains$list

    ## Concatenate total number of infections per MCMC sample
    total_inf_chain <- get_total_number_infections(inf_chain, FALSE)
    chain <- merge(chain, total_inf_chain)

    ## Combine total number of infections with theta chain
  list_total_inf_chains <- lapply(inf_list_chains, function(x) get_total_number_infections(x, FALSE))
  theta_list_chains <- Map(
    function(x, y) merge(data.table(x), data.table(y), by = c("samp_no", "chain_no")),
    theta_list_chains, list_total_inf_chains
  )
  theta_list_chains <- lapply(theta_list_chains, function(x) x[order(x[, "samp_no"]), ])
  unique_samp_nos <- lapply(theta_list_chains, function(x) unique(x[, "samp_no"])$samp_no)
  unique_samp_nos <- Reduce(intersect, unique_samp_nos)
  theta_list_chains <- lapply(theta_list_chains, function(x) x[x$samp_no %in% unique_samp_nos, ])

    ## Convert to MCMC objects if desired
  if (convert_mcmc) {
    theta_list_chains <- lapply(theta_list_chains, as.mcmc)
    chain <- as.mcmc(chain)
  }

  return(list("theta_chain" = chain, "inf_chain" = inf_chain, "theta_list_chains" = theta_list_chains, "inf_list_chains" = inf_list_chains))
}

#' Load MCMC chains for theta
#'
#' Searches the given working directory for MCMC outputs from \code{\link{serosolver}}, loads these in, subsets for burn in and thinning, and formats as both lists and a combined data frame.
#' @param location defaults to current working directory. Gives relative file path to look for files ending in "_chain.csv"
#' @param par_tab if not NULL, can use this to only extract free model parameters
#' @param estimated_only if TRUE, only returns free model parameters (par_tab$fixed == 0) if par_tab specified
#' @param thin thin the chains by every thin'th sample
#' @param burnin discard the first burnin samples from the MCMC chain
#' @param convert_mcmc if TRUE, converts everything to MCMC objects (from the `coda` R package)
#' @return a list with a) a list of each chain separately; b) a combined data frame, indexing each iteration by which chain it comes from
#' @family load_data_functions
#' @examples
#' \dontrun{load_theta_chains(par_tab=par_tab, estimated_only=TRUE,thin=10,burnin=5000,convert_mcmc=TRUE)}
#' @export
load_theta_chains <- function(location = getwd(), par_tab = NULL, estimated_only = TRUE, thin = 1, burnin = 0, convert_mcmc = TRUE) {
  chains <- Sys.glob(file.path(location, "*_chain.csv"))
  message("Chains detected: ", length(chains), sep = "\t")
  if (length(chains) < 1) {
      message("Error - no chains found")
      return(NULL)
  }

  ## Read in the MCMC chains with fread for speed
  read_chains <- lapply(chains, read.csv)

  message("Highest MCMC sample interations: \n")
  lapply(read_chains, function(x) message(max(x$samp_no)))
  
  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[seq(1, nrow(x), by = thin), ])
  read_chains <- lapply(read_chains, function(x) x[x$samp_no > burnin, ])
  max_samp_no <- min(as.numeric(lapply(read_chains, function(x) max(x$samp_no))))
  read_chains <- lapply(read_chains, function(x) x[x$samp_no <= max_samp_no, ])
  unique_samp_nos <- lapply(read_chains, function(x) unique(x[, "samp_no"]))
  unique_samp_nos <- Reduce(intersect, unique_samp_nos)
  read_chains <- lapply(read_chains, function(x) x[x$samp_no %in% unique_samp_nos, ])

  for (i in 1:length(read_chains)) read_chains[[i]]$chain_no <- i

  ## Get the estimated parameters only
  if (estimated_only & !is.null(par_tab)) {
    fixed <- par_tab$fixed
    use_colnames <- intersect(c("samp_no", par_tab$names[which(fixed == 0)], "posterior_prob", "likelihood", "prior_prob", "chain_no"), colnames(read_chains[[1]]))
    read_chains <- lapply(read_chains, function(x) x[, use_colnames])
  }

  ## Try to create an MCMC list. This might not work, which is why we have a try catch
  list_chains <- read_chains
  if (convert_mcmc) {
    list_chains <- tryCatch({
      tmp_list <- lapply(read_chains, coda::as.mcmc)
      tmp_list <- coda::as.mcmc.list(tmp_list)
    },
    warning = function(w) {
      message(w)
      NULL
    }, error = function(e) {
      message(e)
      NULL
    },
    finally = {
      tmp_list
    }
    )
  }

  chain <- do.call("rbind", read_chains)
  if (convert_mcmc) chain <- as.mcmc(chain)
  return(list("list" = list_chains, "chain" = chain))
}

#' Load MCMC chains for infection histories
#'
#' Searches the given working directory for MCMC outputs from \code{\link{serosolver}}, loads these in, subsets for burn in and thinning, and formats as both lists and a combined data table.
#' @param location defaults to current working directory. Where to look for MCMC chains? These are files ending in "_infection_histories.csv"
#' @inheritParams load_theta_chains
#' @param chain_subset if not NULL, a vector of indices to only load and store a subset of the chains detected. eg. chain_subset = 1:3 means that only the first 3 detected files will be processed.
#' @return a list with a) a list of each chain as a data table separately; b) a combined data table, indexing each iteration by which chain it comes from
#' @family load_data_functions
#' @examples
#' \dontrun{load_infection_chains(thin=10,burnin=5000,chain_subset=1:3)}
#' @export
load_infection_chains <- function(location = getwd(), thin = 1, burnin = 0, chain_subset = NULL) {
  chains <- Sys.glob(file.path(location, "*_infection_histories.csv"))
  chains_old <- Sys.glob(file.path(location, "*_infectionHistories.csv"))
  chains <- c(chains, chains_old)
  if (!is.null(chain_subset)) chains <- chains[chain_subset]
  message("Chains detected: ", chains, sep = "\n")
  if (length(chains) < 1) {
    message("Error - no chains found")
    return(NULL)
  }

  message("Reading in infection history chains. May take a while.")
  ## Read in the MCMC chains with fread for speed
  read_chains <- lapply(chains, data.table::fread)

  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[samp_no > burnin, ])
  read_chains <- lapply(read_chains, function(x) {
    samp_nos <- unique(x$samp_no)
    samp_nos <- samp_nos[seq(1, length(samp_nos), by = thin)]
    x <- x[samp_no %in% samp_nos, ]
  })
  max_samp_no <- min(as.numeric(lapply(read_chains, function(x) max(x$samp_no))))
  read_chains <- lapply(read_chains, function(x) x[samp_no <= max_samp_no, ])

  message("Number of rows: ")
  lapply(read_chains, function(x) message(nrow(x)))

  for (i in 1:length(read_chains)) read_chains[[i]]$chain_no <- i
  chain <- do.call("rbind", read_chains)
  return(list("list" = read_chains, "chain" = chain))
}


#' Get total number of infections
#'
#' Finds the total number of infections for each iteration of an MCMC chain
#' @inheritParams plot_infection_history_chains_time
#' @return a data table
#' @examples
#' \dontrun{
#' inf_chain <- load_infection_chains(thin=10,burnin=5000,chain_subset=1:3)
#' n_infs <- get_total_number_infections(inf_chain$chain, pad_chain=FALSE)
#' }
#' @export
get_total_number_infections <- function(inf_chain, pad_chain = TRUE) {
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  n_strain <- max(inf_chain$j)
  data.table::setkey(inf_chain, "samp_no", "chain_no")
  ## For each individual, how many infections did they have in each sample in total?
  n_inf_chain <- inf_chain[, list(total_infections = sum(x)), by = key(inf_chain)]
}

#' Estimate vector mode
#'
#' @param x the vector to be estimated
#' @return the estimated mode of the given vector of values
#' @examples
#' x <- runif(1000)
#' y <- estimate_mode(x)
#' @export
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

#' Formatted quantiles
#'
#' Given a vector of MCMC samples, generates and formats the desired quantile estimates
#' @param x the vector to summarise
#' @param sig_f how many significant figures to print
#' @param qs the vector of quantiles
#' @param as_text if TRUE, formats nicely as text rather than a vector of numbers
#' @return the formatted quantiles
#' @examples
#' data(example_theta_chain)
#' x <- example_theta_chain$boost_long
#' generate_quantiles(x)
#' @export
generate_quantiles <- function(x, sig_f = 3, qs = c(0.025, 0.5, 0.975), as_text = TRUE) {
  res <- signif(quantile(x, qs), sig_f)
  if (as_text) {
    res <- paste(res[2], " (", res[1], "-", res[3], ")", sep = "")
  }
  return(res)
}


#' Generate antibody level credible intervals
#'
#' Generates credible intervals on antibody levels and infection histories from an MCMC chain output.
#' @param chain the full MCMC chain to generate antibody level trajectories from
#' @param infection_histories the MCMC chain for infection histories
#' @param antibody_data the data frame of antibody level data
#' @param individuals the subset of individuals to generate credible intervals for
#' @param antigenic_map (optional) a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. See \code{\link{example_antigenic_map}}
#' @param possible_exposure_times (optional) if no antigenic map is specified, this argument gives the vector of times at which individuals can be infected
#' @param par_tab the table controlling the parameters in the MCMC chain
#' @param nsamp number of samples to take from posterior
#' @param add_residuals if true, returns an extra output summarising residuals between the model prediction and data
#' @param measurement_bias default NULL, optional data frame giving the index of `rho` that each biomarker_id and biomarker_group which uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clusters
#' @param for_res_plot TRUE/FALSE value. If using the output of this for plotting of residuals, returns the actual data points rather than summary statistics
#' @param expand_antibody_data TRUE/FALSE value. If TRUE, solves antibody level predictions for the entire study period (i.e., between the range of antibody_data$sample_time). If left FALSE, then only solves for the infections times at which a antibody level against the circulating biomarker_id was measured in antibody_data.
#' @param expand_to_all_times TRUE/FALSE value. If TRUE, solves antibody level predictions for all possible infection times (i.e., for the range in possible_exposure_times). If left FALSE, then only solves for the infections times at which a antibody level against the circulating biomarker_id was measured in antibody_data.
#' @param antibody_level_before_infection TRUE/FALSE value. If TRUE, solves antibody level predictions, but gives the predicted antibody level at a given time point BEFORE any infection during that time occurs.
#' @param for_regression if TRUE, returns posterior draws rather than posterior summaries
#' @param data_type integer, currently accepting 1 or 2. Set to 1 for discretized, bounded data, or 2 for continuous, bounded data. 
#' @return a list with the antibody level predictions (95% credible intervals, median and multivariate posterior mode) and the probabilities of infection for each individual in each epoch
#' @examples
#' \dontrun{
#' data(example_theta_chain)
#' data(example_inf_chain)
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' data(example_par_tab)
#'
#' y <- get_antibody_level_predictions(example_theta_chain, example_inf_chain, example_antibody_data,
#'                           unique(example_antibody_data$individual), example_antigenic_map,
#'                           example_par_tab,expand_antibody_data = FALSE)
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
                                           antibody_level_before_infection=FALSE, for_regression=FALSE,
                                           data_type=1,start_level="none"){
  par_tab <- add_scale_pars(par_tab,antibody_data,demographics)
  
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
  
  
  if (expand_antibody_data) {
    if(expand_to_all_times){
      expand_times <- possible_exposure_times
    } else {
      expand_times <- unique(antibody_data$sample_time)
    }
    antibody_data1 <- expand.grid(
      individual = unique(antibody_data$individual),
      sample_time = expand_times,
      biomarker_group=unique(antibody_data$biomarker_group),
      measurement = 0, repeat_number = 1
    )
    antibody_data2 <- antibody_data %>% dplyr::select(-c(sample_time,measurement,repeat_number)) %>% distinct()
    antibody_data1 <- merge(antibody_data1, antibody_data2)
    antibody_data1 <- antibody_data1 %>% arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
  }
  model_func <- create_posterior_func(par_tab, antibody_data1, antigenic_map, possible_exposure_times,
                                     prior_version=2,
                                     measurement_bias = measurement_bias, function_type = 4,
                                      antibody_level_before_infection=antibody_level_before_infection,
                                      data_type=data_type,start_level=start_level,
                                     demographics=demographics
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

#' Get posterior information infection histories
#'
#' Finds the median, mean and 95% credible intervals for the attack rates and total number of infections per individual
#' @param solve_cumulative if TRUE, finds the cumulative infection histories for each individual. This takes a while, so is left FALSE by default.
#' @inheritParams plot_posteriors_infhist
#' @return a list of data frames with summary statistics
#' @family infection_history_plots
#' @examples
#' data(example_inf_chain)
#' data(example_antigenic_map)
#' data(example_antibody_data)
#' data(example_inf_hist)
#' possible_exposure_times <- example_antigenic_map$inf_times
#' ## Find number alive in each time period
#' n_alive <- get_n_alive(example_antibody_data, possible_exposure_times)
#' ## Get actual number of infections per time
#' n_infs <- colSums(example_inf_hist)
#' ## Create data frame of true ARs
#' known_ar <- n_infs/n_alive
#' known_ar <- data.frame("j"=possible_exposure_times,"AR"=known_ar,"population_group"=1)
#' ## Get true infection histories
#' known_inf_hist <- data.frame(example_inf_hist)
#' colnames(known_inf_hist) <- possible_exposure_times
#'
#' ## Need to get population_group specific n_alive and adjust to correct time frame 
#' n_alive_group <- get_n_alive_group(example_antibody_data, possible_exposure_times,melt_dat = TRUE)
#' n_alive_group$j <- possible_exposure_times[n_alive_group$j]
#' results <- calculate_infection_history_statistics(example_inf_chain, 0, possible_exposure_times,
#'                                                   n_alive=n_alive_group, known_ar=known_ar,
#'                                                   known_infection_history=known_inf_hist)
#' @export
calculate_infection_history_statistics <- function(inf_chain, burnin = 0, possible_exposure_times = NULL,
                                                   n_alive = NULL, known_ar = NULL,
                                                   group_ids = NULL,
                                                   known_infection_history = NULL,
                                                   solve_cumulative=FALSE) {
  inf_chain <- inf_chain[inf_chain$samp_no > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  ## message("Padding inf chain...\n")
  inf_chain <- pad_inf_chain(inf_chain)
  ## message("Done\n")
  
  if (!is.null(group_ids)) {
    inf_chain <- merge(inf_chain, data.table(group_ids))
  } else {
    inf_chain$population_group <- 1
  }
  
  ## message("Calculating by time summaries...\n")
  data.table::setkey(inf_chain, "population_group", "j", "samp_no", "chain_no")
  n_inf_chain <- inf_chain[, list(total_infs = sum(x)), by = key(inf_chain)]
  
  
  if (!is.null(possible_exposure_times)) {
    n_inf_chain$j <- possible_exposure_times[n_inf_chain$j]
  }
  
  if (!is.null(n_alive)) {
    n_inf_chain <- merge(n_inf_chain, n_alive, by = c("j", "population_group"))
    n_inf_chain$total_infs <- n_inf_chain$total_infs / n_inf_chain$n_alive
    n_inf_chain[is.nan(n_inf_chain$total_infs), "total_infs"] <- 0
  }
  data.table::setkey(n_inf_chain, "population_group", "samp_no", "chain_no")
  n_inf_chain[, cumu_infs := cumsum(total_infs), by = key(n_inf_chain)]
  
  if(length(unique(n_inf_chain$chain_no)) > 1){
    gelman_res_j <- ddply(n_inf_chain, .(population_group,j), function(tmp_chain){
      tmp_chain_mcmc <- split(as.data.table(tmp_chain), by=c("chain_no"))
      tmp_chain_mcmc <- lapply(tmp_chain_mcmc, function(x) as.mcmc(x[,c("total_infs")]))
      tmp_chain_mcmc <- as.mcmc.list(tmp_chain_mcmc)
      gelman.diag(tmp_chain_mcmc)[[1]][1,]
    })
    colnames(gelman_res_j) <- c("population_group","j","gelman_point","gelman_upper")
  } else {
    gelman_res_j <- NULL
  }
  
  
  data.table::setkey(n_inf_chain, "j", "population_group")
  n_inf_chain_summaries <- n_inf_chain[, list(
    mean = mean(as.double(total_infs)), median = median(as.double(total_infs)),
    lower_quantile = quantile(as.double(total_infs), c(0.025)),
    upper_quantile = quantile(as.double(total_infs), c(0.975)),
    effective_size = tryCatch({
      coda::effectiveSize(total_infs)
    }, error = function(e) {
      0
    })
  ),
  by = key(n_inf_chain)
  ]
  if(length(unique(n_inf_chain$chain_no)) > 1){
    
    n_inf_chain_summaries <- merge(n_inf_chain_summaries, gelman_res_j, by=c("j","population_group"))
  }
  n_inf_chain_summaries_cumu <- n_inf_chain[, list(
    mean = mean(as.double(cumu_infs)), median = median(as.double(cumu_infs)),
    lower_quantile = quantile(as.double(cumu_infs), c(0.025)),
    upper_quantile = quantile(as.double(cumu_infs), c(0.975)),
    effective_size = tryCatch({
      coda::effectiveSize(cumu_infs)
    }, error = function(e) {
      0
    })
  ),
  by = key(n_inf_chain)
  ]
  ## message("Done\n")
  if (!is.null(known_ar)) {
    n_inf_chain_summaries <- merge(n_inf_chain_summaries, known_ar, by = c("j","population_group"))
    n_inf_chain_summaries$correct <- (n_inf_chain_summaries$AR >=
                                        n_inf_chain_summaries$lower_quantile) & (n_inf_chain_summaries$AR <=
                                                                                   n_inf_chain_summaries$upper_quantile)
  }
  ## message("Calculating by individual summaries...\n")
  data.table::setkey(inf_chain, "i", "samp_no", "chain_no")
  n_inf_chain_i <- inf_chain[, list(total_infs = sum(x)), by = key(inf_chain)]
  
  data.table::setkey(n_inf_chain_i, "i")
  n_inf_chain_i_summaries <- n_inf_chain_i[, list(
    mean = mean(total_infs),
    median = as.integer(median(total_infs)),
    lower_quantile = quantile(total_infs, 0.025),
    upper_quantile = quantile(total_infs, 0.975),
    effective_size = tryCatch({
      coda::effectiveSize(total_infs)
    }, error = function(e) {
      0
    })
  ),
  by = key(n_inf_chain_i)
  ]
  
  if(solve_cumulative){
    data.table::setkey(inf_chain,"i", "samp_no", "chain_no")
    n_inf_chain_i_cumu <- inf_chain[, cumu_infs := cumsum(x), by = key(inf_chain)]
    
    data.table::setkey(n_inf_chain_i_cumu, "i","j")
    n_inf_chain_i_summaries_cumu <- n_inf_chain_i_cumu[, list(
      mean = mean(cumu_infs),
      median = as.integer(median(cumu_infs)),
      lower_quantile = quantile(cumu_infs, 0.025),
      upper_quantile = quantile(cumu_infs, 0.975),
      effective_size = tryCatch({
        coda::effectiveSize(cumu_infs)
      }, error = function(e) {
        0
      })
    ),
    by = key(n_inf_chain_i_cumu)
    ]
  } else {
    n_inf_chain_i_summaries_cumu <- NULL
  }
 ##  message("Done\n")
  if (!is.null(known_infection_history)) {
    true_n_infs <- rowSums(known_infection_history)
    true_n_infs <- data.frame(i = 1:length(true_n_infs), true_infs = true_n_infs)
    n_inf_chain_i_summaries <- merge(n_inf_chain_i_summaries, true_n_infs, by = "i")
    n_inf_chain_i_summaries$correct <- (n_inf_chain_i_summaries$true_inf >=
                                          n_inf_chain_i_summaries$lower_quantile) & (n_inf_chain_i_summaries$true_inf <=
                                                                                       n_inf_chain_i_summaries$upper_quantile)
  }
  
  return(list(
    "by_year" = n_inf_chain_summaries, "by_indiv" = n_inf_chain_i_summaries,
    "by_year_cumu" = n_inf_chain_summaries_cumu, "by_indiv_cumu" = n_inf_chain_i_summaries_cumu
  ))
}

#' Identify runs of infections from posterior
#'
#' For each individual and MCMC iteration, uses the infection history MCMC chain and detects runs of consecutive infections.
#' @param inf_chain data table of the infection histories posterior
#' @return a tibble giving the consecutive infection run length, the start and end time of each run, which index the run is (ie., which distinct infection), and the time from the end of the previous run, for each i and samp_no
#' @examples
#' \dontrun{
#' identify_run_lengths(inf_chain)
#' }
#' @export
identify_run_lengths <- function(inf_chain) {
  inf_chain %>% 
    arrange(samp_no, i, j) %>%
    as_tibble() %>%
    group_by(i, samp_no) %>%
    mutate(run_group = cumsum(c(0, diff(x != 1) != 0))) %>%
    filter(x == 1) %>%
    group_by(i, samp_no, run_group) %>%
    summarise(run_length = n(),
              start_time = first(j),
              end_time = last(j))%>%
    group_by(i, samp_no) %>%
    mutate(infection_index = row_number(),
           time_diff = start_time - lag(end_time,1,default=NA)) %>%
    ungroup()%>%
    select(-run_group)
}

#' Summarize runs of infections from posterior
#'
#' Either takes the output of \code{\link{identify_run_lengths}} and produces summary statistics giving the posterior median and 95% quantiles for the length of each infection run
#' @param inf_chain data table of the infection histories posterior or tibble from \code{\link{identify_run_lengths}}
#' @return a tibble summarizing the infection run lengths for each individual and distinct infection event
#' @examples
#' \dontrun{
#' summarize_run_lengths(inf_chain)
#' summarize_run_lengths(identify_run_lengths(inf_chain))
#' }
#' @export
summarize_run_lengths <- function(inf_chain){
  if(!("run_length" %in% colnames(inf_chain))){
    summary <- identify_run_lengths(inf_chain) 
  } else {
    summary <- inf_chain
  }
  summary %>%
    group_by(i, infection_index) %>% 
    dplyr::summarize(median_run_length=median(run_length),
                     lower95_run_length=quantile(run_length,0.025),
                     upper95_run_length=quantile(run_length,0.975))
}
