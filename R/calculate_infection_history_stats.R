#' Get posterior information for estimated infection histories
#'
#' Finds the median, mean and 95% credible intervals for the attack rates and total number of infections per individual
#' @param solve_cumulative if TRUE, finds the cumulative infection histories for each individual. This takes a while, so is left FALSE by default.
#' @inheritParams plot_infection_history_posteriors
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
calculate_infection_history_statistics <- function(inf_chain, 
                                                   burnin = 0, 
                                                   possible_exposure_times = NULL,
                                                   n_alive = NULL, known_ar = NULL,
                                                   group_ids = NULL,
                                                   known_infection_history = NULL,
                                                   solve_cumulative=FALSE,
                                                   pad_chain=FALSE) {
  inf_chain <- inf_chain[inf_chain$samp_no > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  ## message("Padding inf chain...\n")
  if(pad_chain){
    inf_chain <- pad_inf_chain(inf_chain)
  }
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
    gelman_res_j <- n_inf_chain %>% 
      dplyr::select(population_group,j,chain_no,total_infs) %>%
      dplyr::group_by(population_group,j,chain_no) %>%
      dplyr::summarize(x = list(as.mcmc(total_infs))) %>%
      dplyr::group_by(population_group,j) %>%
      dplyr::summarize(gelman_point = unlist(gelman.diag(as.mcmc.list(x))[[1]][1,1]),
                       gelman_upper=unlist(gelman.diag(as.mcmc.list(x))[[1]][1,2])) %>%
      ungroup()
    
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
