#' Plot inferred posteriors infection histories
#'
#' Plots and calculates many summary statistics from the infection history MCMC chain
#' @param inf_chain the data table with infection history samples from \code{\link{serosolver}}
#' @param possible_exposure_times vector of the epochs of potential circulation
#' @param n_alive_group vector with the number of people alive in each year of circulation.
#' @param known_ar data frame of known attack rates, if known.
#' @param known_infection_history data frame of known infection histories.
#' @param burnin if not already discarded, discard burn in from chain (takes rows where samp_no > burnin)
#' @param samples how many samples from the chain to take
#' @param pad_chain if TRUE, pads the infection history MCMC chain with non-infection events
#' @return a list of ggplot objects and data frame of posterior estimates
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' ## Load in example data
#' data(example_inf_chain)
#' data(example_antigenic_map)
#' data(example_antibody_data)
#'
#' possible_exposure_times <- example_antigenic_map$inf_times
#' ## Setup known attack rates
#' n_alive <- get_n_alive(example_antibody_data, possible_exposure_times)
#' n_infs <- colSums(example_inf_hist)
#' known_ar <- n_infs/n_alive
#' known_ar <- data.frame("j"=possible_exposure_times,"AR"=known_ar,"population_group"=1)
#' 
#' ## Setup known infection histories
#' known_inf_hist <- data.frame(example_inf_hist)
#' colnames(known_inf_hist) <- possible_exposure_times
#' 
#' n_alive_group <- get_n_alive_group(example_antibody_data, possible_exposure_times,melt_dat = TRUE)
#' n_alive_group$j <- possible_exposure_times[n_alive_group$j]
#' all_plots <- plot_infection_history_posteriors(example_inf_chain, possible_exposure_times, n_alive_group,
#'                                      known_ar=known_ar,known_infection_history = known_inf_hist,
#'                                      samples=100)
#' }
#' @export
plot_infection_history_posteriors <- function(inf_chain,
                                    possible_exposure_times,
                                    n_alive,
                                    known_ar = NULL,
                                    known_infection_history = NULL,
                                    burnin = 0,
                                    samples = 100,
                                    pad_chain = FALSE) {
  ## Discard burn in period if necessary
  inf_chain <- inf_chain[inf_chain$samp_no > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if(samples > length(unique(inf_chain$samp_no))) {
    samples <- min(samples, length(unique(inf_chain$samp_no)))
    message("Number of samples requested is greater than number of MCMC samples provided, reducing to maximum number of possible samples")
  }
  
  ## Pads the chain with 0 entries (ie. non-infection events)
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  
  ## Thin the chain a bit to increase speed for plots
  thin_chain <- inf_chain[inf_chain$samp_no %in% sample(unique(inf_chain$samp_no), samples, replace = FALSE), ]
  
  ## Trace plot by time
  n_time_sample <-  min(length(possible_exposure_times),10)
  time_plot <- plot_infection_history_chains_time(thin_chain, 0, sample(1:length(possible_exposure_times), n_time_sample), n_alive=NULL, FALSE)
  ## Trace plot by indiv
  n_indiv_sample <- min(length(unique(inf_chain$i)),10)
  indiv_plot <- plot_infection_history_chains_indiv(thin_chain, 0, sample(unique(inf_chain$i),n_indiv_sample,replace=FALSE), FALSE)
  ## Pointrange plot for total number of infections
  number_plot <- plot_individual_number_infections(thin_chain, FALSE)
  
  ## Posterior summaries and ESS
  results <- calculate_infection_history_statistics(inf_chain, 0, possible_exposure_times,
                                                    n_alive, 
                                                    known_ar=known_ar,
                                                    known_infection_history=known_infection_history,
                                                    pad_chain=pad_chain
  )
  return(list(
    "by_time_trace" = time_plot, "by_indiv_trace" = indiv_plot,
    "indiv_infections" = number_plot, "estimates" = results
  ))
}




#' Plot MCMC trace for infections per year
#'
#' @param inf_chain the data table with infection history samples from \code{\link{serosolver}}
#' @param burnin optionally remove all samp_no < burnin from the chain
#' @param times vector of integers, if not NULL, only plots a subset of years (where 1 is the first year eg. 1968)
#' @param n_alive if not NULL, then divides number of infections per year by number alive to give attack rates rather than total infections
#' @param pad_chain if TRUE, pads the infection history MCMC chain to have entries for non-infection events
#' @return a list of two ggplot objects - the MCMC trace and MCMC densities
#' @seealso \code{\link{plot_infection_history_chains_indiv}}
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_inf_chain)
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_times
#' n_alive_group <- get_n_alive_group(example_antibody_data, possible_exposure_times,melt_dat = TRUE)
#' n_alive_group$j <- possible_exposure_times[n_alive_group$j]
#' plot_infection_history_chains_time(example_inf_chain, 0, sample(1:length(times),10),n_alive,FALSE)
#' }
#' @export
plot_infection_history_chains_time <- function(inf_chain, burnin = 0, times = NULL,
                                               n_alive = NULL, pad_chain = TRUE) {
  inf_chain <- inf_chain[inf_chain$samp_no > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  data.table::setkey(inf_chain, "j", "samp_no", "chain_no")
  n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
  
  if (!is.null(n_alive)) {
    n_inf_chain$V1 <- n_inf_chain$V1 / n_alive[match(n_inf_chain$j, n_alive$j),"n_alive"]
  }
  
  if (!is.null(times)) {
    use_times<- intersect(unique(n_inf_chain$j), times)
    n_inf_chain <- n_inf_chain[n_inf_chain$j %in% use_times, ]
  }
  n_inf_chain$chain_no <- as.factor(n_inf_chain$chain_no)
  inf_chain_p <- ggplot(n_inf_chain) + geom_line(aes(x = samp_no, y = V1, col = chain_no)) +
    ylab("Estimated attack rate") +
    xlab("MCMC sample") +
    theme_bw() +
    facet_wrap(~j, scales = "free_y")
  inf_chain_den <- ggplot(n_inf_chain) + geom_density(aes(x = V1, fill = chain_no),alpha=0.5) +
    xlab("Estimated attack rate") +
    ylab("Posterior density") +
    theme_bw() +
    facet_wrap(~j, scales = "free_x")
  return(list(inf_chain_p, inf_chain_den))
}

#' Plot MCMC trace for infections per individual
#'
#' @inheritParams plot_infection_history_chains_time
#' @param indivs vector of integers, if not NULL, only plots a subset of individuals (where 1 is the first individual)
#' @return a list of two ggplot objects - the MCMC trace and MCMC densities
#' @seealso \code{\link{plot_infection_history_chains_indiv}}
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_inf_chain)
#' plot_infection_history_chains_indiv(example_inf_chain, 0, 1:10, FALSE)
#' }
#' @export
plot_infection_history_chains_indiv <- function(inf_chain, burnin = 0, indivs = NULL, pad_chain = TRUE) {
  inf_chain <- inf_chain[inf_chain$samp_no > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  data.table::setkey(inf_chain, "i", "samp_no", "chain_no")
  n_inf_chain_i <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
  
  if (!is.null(indivs)) {
    use_indivs <- intersect(unique(n_inf_chain_i$i), indivs)
    n_inf_chain_i <- n_inf_chain_i[n_inf_chain_i$i %in% use_indivs, ]
  }
  n_inf_chain_i$chain_no <- as.factor(n_inf_chain_i$chain_no)
  inf_chain_p_i <- ggplot(n_inf_chain_i) + geom_line(aes(x = samp_no, y = V1, col = chain_no)) +
    ylab("Estimated total number of infections") +
    xlab("MCMC sample") +
    theme_bw() +
    facet_wrap(~i)
  inf_chain_den_i <- ggplot(n_inf_chain_i) + geom_histogram(aes(x = V1, fill = chain_no), binwidth = 1) +
    xlab("Estimated total number of infections") +
    ylab("Posterior density") +
    theme_bw() +
    facet_wrap(~i)
  return(list(inf_chain_p_i, inf_chain_den_i))
}


#' Total number of infections
#'
#' Plots the total number of inferred infections in the MCMC chain as a trace plot and density plot
#' @inheritParams plot_infection_history_chains_time
#' @return two ggplot2 objects
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_inf_chain)
#' plot_total_number_infections(example_inf_chain)
#' }
#' @export
plot_total_number_infections <- function(inf_chain, pad_chain = TRUE) {
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  n_inf_chain <- get_total_number_infections(inf_chain, pad_chain)
  n_inf_chain$chain_no <- as.factor(n_inf_chain$chain_no)
  inf_chain_p <- ggplot(n_inf_chain) + geom_line(aes(x = samp_no, y = total_infections, col = chain_no)) +
    ylab("Estimated total number of infections") +
    xlab("MCMC sample") +
    theme_bw()
  inf_chain_den <- ggplot(n_inf_chain) + geom_density(aes(x = total_infections, fill = chain_no),alpha=0.5) +
    xlab("Estimated total number of infections") +
    ylab("Posterior density") +
    theme_bw()
  return(list(inf_chain_p, inf_chain_den))
}


#' Plot point range number infections per individual
#'
#' @inheritParams plot_infection_history_chains_time
#' @return a ggplot object
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_inf_chain)
#' plot_individual_number_infections(example_inf_chain)
#' }
#' @export
plot_individual_number_infections <- function(inf_chain, pad_chain = TRUE) {
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  n_strain <- max(inf_chain$j)
  data.table::setkey(inf_chain, "i", "samp_no", "chain_no")
  ## For each individual, how many infections did they have in each sample in total?
  n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
  
  ## Get quantiles on total number of infections per indiv across all samples
  indiv_hist <- n_inf_chain %>% group_by(i) %>% 
    dplyr::summarize(lower=quantile(V1,0.025),median=quantile(V1,0.5),upper=quantile(V1,0.975)) %>% ungroup() %>% dplyr::rename(individual=i)
  indiv_hist <- indiv_hist[order(indiv_hist$median), ]
  indiv_hist$individual <- 1:nrow(indiv_hist)
  p <- ggplot(indiv_hist) +
    geom_pointrange(aes(x = individual + 1, y = median, ymin = lower, ymax = upper),
                    linewidth = 0.1, shape = 21, fatten = 0.1
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Individual (ordered)") +
    ylab("Estimated number of infections") +
    theme_bw()
  return(p)
}

#' Plot cumulative and per time posterior infection probability densities
#'
#' For each individual requested, plots the median and 95% quantiles on a) the cumulative number of infections over a lifetime and b) the posterior probability that an infection occured in a given time point
#' @param inf_chain the infection history chain
#' @param burnin only plot samples where samp_no > burnin
#' @param indivs vector of individual ids to plot
#' @param real_inf_hist if not NULL, adds lines to the plots showing the known true infection times
#' @param start_inf if not NULL, adds lines to show where the MCMC chain started
#' @param possible_exposure_times vector of times at which individuals could have been infected
#' @param nsamp how many samples from the MCMC chain to take?
#' @param ages if not NULL, adds lines to show when an individual was born
#' @param number_col how many columns to use for the cumulative infection history plot
#' @param pad_chain if TRUE, pads the infection history MCMC chain to have entries for non-infection events
#' @param subset_times if not NULL, pass a vector of indices to only take a subset of indices from possible_exposure_times
#' @param return_data if TRUE, returns the infection history posterior densities used to generate the plots
#' @return two ggplot objects
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_inf_chain)
#' data(example_antigenic_map)
#' data(example_inf_hist)
#' data(example_antibody_data)
#' 
#' ages <- unique(example_antibody_data[,c("individual","birth")])
#' times <- example_antigenic_map$inf_times
#' indivs <- 1:10
#' plot_cumulative_infection_histories(example_inf_chain, 0, indivs, example_inf_hist, NULL, times,
#'                               ages=ages, number_col=2,pad_chain=FALSE, return_data=TRUE)
#' }
#' @export
plot_cumulative_infection_histories <- function(inf_chain, burnin = 0, indivs, real_inf_hist = NULL, start_inf = NULL,
                                                possible_exposure_times, nsamp = 100, ages = NULL, number_col = 1,
                                                pad_chain = TRUE, subset_times = NULL, return_data = FALSE) {
  inf_chain <- inf_chain[inf_chain$samp_no > burnin, ]
  inf_chain <- inf_chain[inf_chain$i %in% indivs,]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain,times=seq_along(possible_exposure_times),indivs=indivs)
  
  ## Create samp_no/chain_no combo and use for samples
  samps_all <- inf_chain[,c("samp_no","chain_no")] %>% distinct() %>% mutate(n = 1:n())
  samps <- sample(unique(samps_all$n),min(nsamp,nrow(samps_all)))
  inf_chain <- inf_chain %>% left_join(samps_all,by=c("samp_no","chain_no")) %>% filter(n %in% samps) %>% select(-n)
  #samps <- sample(unique(inf_chain$samp_no), nsamp)
  #inf_chain <- inf_chain[inf_chain$samp_no %in% samps, ]
  
  ## Get number of probability that infected in a given time point per individual and year
  inf_chain1 <- inf_chain[inf_chain$i %in% indivs, ]
  inf_chain1$j <- as.numeric(possible_exposure_times[inf_chain1$j])
  if (!is.null(subset_times)) inf_chain1 <- inf_chain1[inf_chain1$j %in% subset_times, ]
  data.table::setkey(inf_chain1, "i", "j")
  # max_samp_no <- length(unique(inf_chain1$samp_no))
  max_samp_no <- nsamp
  
  ## Number of samples with a 1 divided by total samples
  densities <- inf_chain1[, list(V1 = sum(x) / max_samp_no), by = key(inf_chain1)]
  
  ## Convert to real time points
  densities$i <- as.numeric(densities$i)
  
  ## If someone wasn't infected in a given year at all, then need a 0
  possible_exposure_times1 <- possible_exposure_times
  if (!is.null(subset_times)) possible_exposure_times1 <- possible_exposure_times[subset_times]
  all_combos <- data.table(expand.grid(i = indivs, j = possible_exposure_times1))
  all_combos$j <- as.numeric(all_combos$j)
  all_combos$i <- as.numeric(all_combos$i)
  all_combos <- data.table::fsetdiff(all_combos[, c("i", "j")], densities[, c("i", "j")])
  all_combos$V1 <- 0
  densities <- rbind(all_combos, densities)
  infection_history1 <- NULL
  if (!is.null(real_inf_hist)) {
    infection_history1 <- as.data.frame(real_inf_hist)
    infection_history1 <- infection_history1[indivs, ]
    infection_history1$individual <- indivs
    
    colnames(infection_history1) <- c(possible_exposure_times, "i")
    infection_history1 <- reshape2::melt(infection_history1, id.vars = "i")
    infection_history1$variable <- as.numeric(as.character(infection_history1$variable))
    infection_history1 <- infection_history1[infection_history1$value == 1, ]
    if (!is.null(subset_times)) infection_history1 <- infection_history1[infection_history1$variable %in% subset_times,]
    
  }
  densities$chain_no <- as.factor(densities$chain_no)
  density_plot <- ggplot() +
    geom_ribbon(data = densities, aes(x = j, ymin=0,ymax = V1),fill="orange")
  if (!is.null(real_inf_hist)) {
    density_plot <- density_plot +
      geom_vline(data = infection_history1, aes(xintercept = variable), col = "black",linetype="dashed",size=0.75)
  }
  density_plot <- density_plot +
    facet_wrap(~i) +
    xlab("Time") +
    ylab("Density") +
    theme(legend.position = "bottom") +
    theme_bw()
  
  ## Generate lower, upper and median cumulative infection histories from the
  ## MCMC chain
  
  tmp_inf_chain <- inf_chain1[inf_chain1$i %in% indivs, ]
  hist_profiles <- tmp_inf_chain %>% arrange(i,samp_no,chain_no,j) %>% group_by(i, samp_no,chain_no) %>% mutate(cumu_x=cumsum(x))
  quant_hist <- hist_profiles %>% group_by(i,j) %>% summarize(lower=quantile(cumu_x, 0.025),upper=quantile(cumu_x,0.975),median=quantile(cumu_x,0.5)) %>% ungroup() %>% rename(individual=i,variable=j)
  #quant_hist$variable <- possible_exposure_times[quant_hist$variable]
  ## If available, process the real infection history matrix for plotting
  real_hist_profiles <- NULL
  if (!is.null(real_inf_hist)) {
    real_hist_profiles <- as.data.frame(t(apply(real_inf_hist, 1, cumsum)))
    colnames(real_hist_profiles) <- possible_exposure_times
    
    real_hist_profiles <- real_hist_profiles[indivs, ]
    real_hist_profiles$individual <- indivs
    real_hist_profiles <- reshape2::melt(real_hist_profiles, id.vars = "individual")
    if (!is.null(subset_times)) real_hist_profiles <- real_hist_profiles[real_hist_profiles$variable %in% subset_times,]
    
  }
  
  ## Process starting point from MCMC chain
  if (!is.null(start_inf)) {
    start_hist_profiles <- as.data.frame(t(apply(start_inf, 1, cumsum)))
    colnames(start_hist_profiles) <- possible_exposure_times
    start_hist_profiles <- start_hist_profiles[indivs, ]
    start_hist_profiles$individual <- indivs
    start_hist_profiles <- reshape2::melt(start_hist_profiles, id.vars = "individual")
  }
  #quant_hist$chain_no <- as.factor(quant_hist$chain_no)
  p1 <- ggplot(quant_hist[quant_hist$individual %in% indivs, ]) +
    geom_line(aes(x = as.integer(as.character(variable)), y = median),col="black") +
    geom_ribbon(aes(x = as.integer(as.character(variable)), ymin = lower, ymax = upper),fill="orange", alpha = 0.5)
  
  if (!is.null(real_inf_hist)) {
    p1 <- p1 +
      geom_line(data = real_hist_profiles[real_hist_profiles$individual %in% indivs, ], aes(x = as.integer(as.character(variable)), y = value), col = "black",linetype="dashed",size=0.75)
  }
  if (!is.null(start_inf)) {
    p1 <- p1 +
      geom_line(data = start_hist_profiles[start_hist_profiles$individual %in% indivs, ], aes(x = as.integer(as.character(variable)), y = value), col = "red")
  }
  
  if (!is.null(ages)) {
    tmp_age <- ages[ages$individual %in% indivs, ]
    age_mask <- create_age_mask(tmp_age[, 2], possible_exposure_times)
    age_dat <- data.frame(j = age_mask, individual = indivs[order(indivs)])
    p1 <- p1 + geom_vline(data = age_dat, aes(xintercept = possible_exposure_times[j]), col = "purple", linetype = "dashed")
  }
  p1 <- p1 + facet_wrap(~individual, ncol = number_col) +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylab("Cumulative infections") +
    xlab("Circulation time")
  return_list <- NULL
  if (return_data) {
    return_list <- list(
      "cumu_infections" = p1, "density_plot" = density_plot,
      "density_data" = densities, "real_hist" = real_hist_profiles,
      "cumu_data" = quant_hist
    )
  } else {
    return_list <- list("cumu_infections" = p1, "density_plot" = density_plot)
  }
  return(return_list)
}
