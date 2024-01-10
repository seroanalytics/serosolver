#' Plot inferred posteriors infection histories
#'
#' Plots and calculates many summary statistics from the infection history MCMC chain
#' @param inf_chain the data table with infection history samples from \code{\link{run_MCMC}}
#' @param years vector of the epochs of potential circulation
#' @param n_alive_group vector with the number of people alive in each year of circulation.
#' @param known_ar data frame of known attack rates, if known.
#' @param known_infection_history data frame of known infection histories.
#' @param burnin if not already discarded, discard burn in from chain (takes rows where sampno > burnin)
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
#' all_plots <- plot_posteriors_infhist(example_inf_chain, possible_exposure_times, n_alive_group,
#'                                      known_ar=known_ar,known_infection_history = known_inf_hist,
#'                                      samples=100)
#' }
#' @export
plot_posteriors_infhist <- function(inf_chain,
                                    years,
                                    n_alive,
                                    known_ar = NULL,
                                    known_infection_history = NULL,
                                    burnin = 0,
                                    samples = 100,
                                    pad_chain = TRUE) {
  ## Discard burn in period if necessary
  inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (is.null(inf_chain$population_group)) {
    inf_chain$population_group <- 1
  }
  if(samples > length(unique(inf_chain$sampno))) {
    samples <- min(samples, length(unique(inf_chain$sampno)))
    message("Number of samples requested is greater than number of MCMC samples provided, reducing to maximum number of possible samples")
  }
  
  ## Pads the chain with 0 entries (ie. non-infection events)
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  
  ## Thin the chain a bit to increase speed for plots
  thin_chain <- inf_chain[inf_chain$sampno %in% sample(unique(inf_chain$sampno), samples, replace = FALSE), ]
  
  ## Trace plot by time
  time_plot <- plot_infection_history_chains_time(thin_chain, 0, sample(1:length(years), 10), n_alive=NULL, FALSE)
  ## Trace plot by indiv
  indiv_plot <- plot_infection_history_chains_indiv(thin_chain, 0, sample(unique(inf_chain$i),10,replace=FALSE), FALSE)
  ## Pointrange plot for total number of infections
  number_plot <- plot_number_infections(thin_chain, FALSE)
  
  ## Posterior summaries and ESS
  results <- calculate_infection_history_statistics(inf_chain, 0, years,
                                                    n_alive, 
                                                    known_ar=known_ar,
                                                    known_infection_history=known_infection_history
  )
  return(list(
    "by_time_trace" = time_plot, "by_indiv_trace" = indiv_plot,
    "indiv_infections" = number_plot, "estimates" = results
  ))
}


#' Plot historical attack rates
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories. 
#' @inheritParams plot_attack_rates
#' @param ymax Numeric. the maximum y value to put on the axis. Default = 1.
#' @param resolution Integer. How many buckets of time is each year split into? ie. 12 for monthly data, 4 for quarterly etc. Default = 1 for annual.
#' @param cumulative if TRUE, plots the cumulative attack rate
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates <- function(infection_histories, antibody_data, possible_exposure_times,
                                      n_alive = NULL, ymax = 1, buckets = 1,
                                      pad_chain = TRUE, true_ar = NULL, by_group = FALSE, group_subset = NULL,
                                      cumulative = FALSE,add_box=FALSE) {
  if (is.null(infection_histories$chain_no)) {
    infection_histories$chain_no <- 1
  }
  if (is.null(infection_histories$population_group)) {
    infection_histories$population_group <- 1
  }
  ## Some year/sample combinations might have no infections there.
  ## Need to make sure that these get considered
  if (pad_chain) infection_histories <- pad_inf_chain(infection_histories)
  ## Subset of groups to plot
  if (is.null(group_subset)) {
    group_subset <- unique(antibody_data$population_group)
  }
  if (!by_group) antibody_data$population_group <- 1
  
  ## Find inferred total number of infections from the MCMC output
  ## Scale by number of individuals that were alive in each epoch
  ## and generate quantiles
  months <- 1:resolution
  years <- possible_exposure_times
  years <- range(floor(years / resolution))
  years <- years[1]:years[2]
  labels <- c(sapply(years, function(x) paste0(months, "/", x)))
  labels1 <- labels[1:length(possible_exposure_times)]
  labels1 <- labels1[seq(1, length(labels1), by = 1)]
  year_break <- possible_exposure_times[seq(1, length(possible_exposure_times), by = 1)]
  
  if (is.null(n_alive)) {
    n_alive <- as.data.frame(get_n_alive_group(antibody_data, possible_exposure_times))
    n_alive$population_group <- 1:nrow(n_alive)
  }
  n_alive_tmp <- reshape2::melt(n_alive, id.vars = "population_group")
  n_alive_tmp$variable <- as.numeric(n_alive_tmp$variable)
  colnames(n_alive_tmp) <- c("population_group", "j", "n_alive")
  
  colnames(infection_histories)[1] <- "individual"
  infection_histories <- merge(infection_histories, data.table(unique(antibody_data[, c("individual", "population_group")])), by = c("individual","population_group"))
  ## Sum infections per year for each MCMC sample
  data.table::setkey(infection_histories, "sampno", "j", "population_group", "chain_no")
  tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
  
  if (cumulative & !pad_chain) message("Error - cannot calculate cumulative incidence without pad_chain = TRUE\n")
  
  tmp <- merge(tmp, data.table(n_alive_tmp), by = c("population_group", "j"))
  tmp$time <- possible_exposure_times[tmp$j]
  tmp$V1 <- tmp$V1 / tmp$n_alive
  tmp[is.nan(tmp$V1), "V1"] <- 0
  if (cumulative && pad_chain) {
    data.table::setkey(tmp, "sampno", "population_group", "chain_no")
    tmp[, V1 := cumsum(V1), by = key(tmp)]
  }
  quantiles <- ddply(tmp, .(time, population_group), function(x) quantile(x$V1, c(0.025,0.25, 0.5,0.75, 0.975)))
  colnames(quantiles) <- c("time", "population_group", "lower","lower2", "median", "upper2","upper")
  
  
  
  p <- ggplot(quantiles[quantiles$population_group %in% group_subset, ])
  if(add_box){
    x_box_min <- min(antibody_data$sample_time)
    p <- p +geom_rect(xmin=x_box_min, xmax=max(quantiles$time),ymin=-1,ymax=2,fill="gray90",alpha=0.1)
  }
  p <- p +
    geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
    geom_ribbon(aes(x = time, ymin = lower2, ymax = upper2), fill = "red", alpha = 0.4) +
    geom_line(aes(x = time, y = median), col = "red")
  if (!is.null(true_ar)) {
    p <- p +
      geom_line(
        data = true_ar[true_ar$population_group %in% group_subset, ], aes(x = time, y = AR),
        col = "purple", size = 0.5
      )
  }
  p <- p +
    ## geom_point(aes(x = year, y = median), col = "purple", size = 0.5) +
    facet_wrap(~population_group, ncol = 2) +
    scale_y_continuous(expand = c(0, 0)) +
    # #scale_x_continuous(expand = c(0, 0), breaks = year_break, labels = labels1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Estimated per capita incidence") +
    xlab("Date")
  return(p)
}

#' Plot historical attack rates with pointrange plots
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories, with pointrange plots for per-time incidence estimates
#' @param infection_histories the MCMC chain for infection histories
#' @param antibody_data the data frame of antibody data
#' @param possible_exposure_times vector of the epochs of potential infection
#' @param n_alive vector with the number of people alive in each year of possible infection Can be left as NULL, and the `birth` variable in `antibody_data` will be used to calculate the number alive
#' @param resolution divides `possible_exposure_times` by this number for x axis labels
#' @param pointsize graphics option, numeric - how big should each point be?
#' @param fatten graphics option, numeric - fatten parameter for ggplot pointrange
#' @param pad_chain if TRUE, fills the infection history data table with entries for non-infection events (ie. 0s). Can be switched to FALSE for speed to get a rough idea of what the attack rates look like.
#' @param prior_pars if not NULL, a list of parameters for the attack rate prior, giving the assumed prior_version along with infection_model_prior_shape1 and infection_model_prior_shape2
#' @param plot_den if TRUE, produces a violin plot of attack rates rather than pointrange
#' @param true_ar data frame of true attack rates, with first column `time` matching `possible_exposure_times`, and second column `AR` giving the attack rate. Column names: population_group, time, AR
#' @param by_group if TRUE, facets the plot by population_group ID
#' @param group_subset if not NULL, plots only this subset of groups eg. 1:5
#' @param plot_residuals if TRUE, plots the residuals between inferred and true attack rate
#' @param colour_by_taken if TRUE, then colours the attack rates by whether or not titres against the circulating antigen at that time were measured
#' @param by_val frequency of x-axis labels
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates_pointrange <- function(infection_histories, antibody_data, possible_exposure_times, 
                              n_alive = NULL,
                              pointsize = 1, fatten = 1,
                              pad_chain = TRUE, prior_pars = NULL,
                              plot_den = FALSE,
                              true_ar = NULL, by_group = FALSE,
                              group_subset = NULL, plot_residuals = FALSE,
                              colour_by_taken = TRUE, by_val = 5,
                              min_time=min(possible_exposure_times),max_time=max(possible_exposure_times)) {
  ## Some year/sample combinations might have no infections there.
  ## Need to make sure that these get considered
  if (is.null(infection_histories$chain_no)) {
    infection_histories$chain_no <- 1
  }
  
  if (pad_chain) infection_histories <- pad_inf_chain(infection_histories)
  
  ## Subset of groups to plot
  if (is.null(group_subset)) {
    group_subset <- unique(antibody_data$population_group)
  }
  if (!by_group) {
    antibody_data$population_group <- 1
    infection_histories$population_group <- 1
  }
  
  ## Find inferred total number of infections from the MCMC output
  ## Scale by number of individuals that were alive in each epoch
  ## and generate quantiles
  if (is.null(n_alive)) {
    n_alive <- get_n_alive_group(antibody_data, possible_exposure_times)
  }
  n_alive <- as.data.frame(n_alive)
  n_alive$population_group <- 1:nrow(n_alive)
  
  n_groups <- length(unique(antibody_data$population_group))
  n_alive_tot <- get_n_alive(antibody_data, possible_exposure_times)
  colnames(infection_histories)[1] <- "individual"
  infection_histories <- merge(infection_histories, data.table(unique(antibody_data[, c("individual", "population_group")])), by = c("individual","population_group"))
  years <- c(possible_exposure_times, max(possible_exposure_times) + 2)
  data.table::setkey(infection_histories, "sampno", "j", "chain_no", "population_group")
  tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
  tmp$taken <- years[tmp$j] %in% unique(antibody_data$sample_time)
  tmp$taken <- ifelse(tmp$taken, "Yes", "No")
  prior_dens <- NULL
  n_alive1 <- n_alive
  if (!is.null(prior_pars)) {
    n_alive$Prior <- 1
    prior_ver <- prior_pars[["prior_version"]]
    alpha1 <- prior_pars[["infection_model_prior_shape1"]]
    beta1 <- prior_pars[["infection_model_prior_shape2"]]
    if (prior_ver == 3) {
      prior_dens <- rbinom(10000, size = max(n_alive_tot), p = alpha1 / (alpha1 + beta1)) / max(n_alive_tot)
    } else {
      prior_dens <- rbeta(10000, alpha1, beta1)
    }
    prior_dens <- data.frame(
      sampno = 1:length(prior_dens), j = max(tmp$j) + 1,
      chain_no = 1, V1 = prior_dens, taken = "Prior", population_group = 1
    )
    prior_dens_all <- NULL
    for (i in 1:nrow(n_alive)) {
      prior_dens$population_group <- i
      prior_dens_all <- rbind(prior_dens_all, prior_dens)
    }
    tmp <- rbind(tmp, prior_dens_all)
  }
  n_alive_tmp <- reshape2::melt(n_alive, id.vars = "population_group")
  n_alive_tmp$variable <- as.numeric(n_alive_tmp$variable)
  colnames(n_alive_tmp) <- c("population_group", "j", "n_alive")
  tmp <- merge(tmp, data.table(n_alive_tmp), by = c("population_group", "j"))
  tmp$V1 <- tmp$V1 / tmp$n_alive
  
  year_breaks <- c(min_time, seq(by_val * round(min_time / by_val), max_time, by = by_val))
  year_labels <- c(min_time, seq(by_val * round(min_time / by_val), max_time, by = by_val))
  
  if (!is.null(prior_dens)) {
    year_breaks <- c(year_breaks, max_time + 2)
    year_labels <- c(year_labels, "Prior")
  }
  
  if (!plot_den) {
    quantiles <- ddply(tmp, .(j, population_group), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
    colnames(quantiles) <- c("j", "population_group", "lower", "median", "upper")
    # quantiles[c("lower", "median", "upper")] <- quantiles[c("lower", "median", "upper")]# / n_alive1[quantiles$j]
    quantiles$j <- years[quantiles$j]
    quantiles$taken <- quantiles$j %in% unique(antibody_data$sample_time)
    quantiles$taken <- ifelse(quantiles$taken, "Yes", "No")
    
    quantiles$tested <- quantiles$j %in% unique(antibody_data$biomarker_id)
    quantiles$tested <- ifelse(quantiles$tested, "Yes", "No")
    
    if (!is.null(prior_dens)) {
      quantiles[quantiles$j == max(years), "taken"] <- "Prior"
    }
    
    colnames(quantiles)[which(colnames(quantiles) == "taken")] <- "Sample taken"
    colnames(quantiles)[which(colnames(quantiles) == "tested")]  <- "Biomarker tested"
    
    p <- ggplot(quantiles[quantiles$population_group %in% group_subset, ]) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(breaks = year_breaks, labels = year_labels) +
      coord_cartesian(ylim=c(0,1)) +
      theme_classic() +
      ylab("Estimated attack rate") +
      xlab("Year")
    
    ## Colour depending on whether or not titres were taken in each year
    if (colour_by_taken == TRUE) {
      p <- p + geom_pointrange(aes(
        x = j, y = median, ymin = lower, ymax = upper,
        col = `Sample taken`, shape = `Sample taken`
      ),
      size = pointsize,
      fatten = fatten
      )
    } else {
      p <- p + geom_pointrange(aes(
        x = j, y = median, ymin = lower, ymax = upper,
        col = `Biomarker tested`, shape = `Biomarker tested`
      ),
      size = pointsize,
      fatten = fatten
      )
    }
  } else {
    tmp$j <- years[tmp$j]
    colnames(tmp)[which(colnames(tmp) == "j"] <- "time"
    if (plot_residuals) {
      true_ar <- true_ar[, c("population_group", "time", "AR")]
      if (!is.null(prior_pars)) {
        true_ar <- rbind(true_ar, data.frame(
          population_group = 1:n_groups,
          time = max(possible_exposure_times) + 3,
          AR = median(prior_dens$V1)
        ))
      }
      tmp <- merge(tmp, true_ar, by = c("population_group", "time"))
      tmp$V1 <- tmp$V1 - tmp$AR
    }
    
    p <- ggplot(tmp[tmp$population_group %in% group_subset, ]) +
      geom_violin(aes(x = time, y = V1, fill = taken, group = j),
                  alpha=0.25,
                  draw_quantiles = c(0.025,0.5,0.975), scale = "width",
                  adjust=2
      )
  }
  if (!is.null(true_ar) & !plot_residuals) {
    p <- p +
      geom_point(
        data = true_ar[true_ar$population_group %in% group_subset, ], aes(x = time, y = AR,shape="True attack rate"),stroke=1.25,
        col = "red", size = 2.5
      )
  }
  if (!plot_residuals) {
    p <- p +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim=c(0,1))
  }
  
  if (by_group) {
    p <- p + facet_wrap(~population_group, ncol = 2)
  }
  if (!is.null(prior_dens)) {
    max_time <- max_time + 2.5
  }
  
  p <- p +
    scale_x_continuous(breaks = year_breaks, labels = year_labels,limits=c(min_time-0.5,max_time)) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ylab("Estimated attack rate") +
    xlab("Date") +
    scale_fill_manual(name="Samples taken in\n time period",values=c("No"="darkorange","Yes"="blue","Prior"="grey40")) +
    scale_shape_manual(name="",values=c("True attack rate"=4))
  
  if (plot_residuals) {
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))
    p_res <- ggplot(tmp) + geom_density(aes(x = V1), fill = "grey40") +
      facet_wrap(~population_group) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "blue") +
      scale_x_continuous(limits = c(-1, 1)) +
      theme_classic() + xlab("Estimated AR - true AR") + ylab("Density")
    return(list(p, p_res))
  }
  
  return(p)
}


#' Plot MCMC trace for infections per year
#'
#' @param inf_chain the data table with infection history samples from \code{\link{run_MCMC}}
#' @param burnin optionally remove all sampno < burnin from the chain
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
  inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  data.table::setkey(inf_chain, "j", "sampno", "chain_no")
  n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
  
  if (!is.null(n_alive)) {
    n_inf_chain$V1 <- n_inf_chain$V1 / n_alive[match(n_inf_chain$j, n_alive$j),"n_alive"]
  }
  
  if (!is.null(years)) {
    use_years <- intersect(unique(n_inf_chain$j), years)
    n_inf_chain <- n_inf_chain[n_inf_chain$j %in% use_years, ]
  }
  
  inf_chain_p <- ggplot(n_inf_chain) + geom_line(aes(x = sampno, y = V1, col = as.factor(chain_no))) +
    ylab("Estimated attack rate") +
    xlab("MCMC sample") +
    theme_bw() +
    facet_wrap(~j, scales = "free_y")
  inf_chain_den <- ggplot(n_inf_chain) + geom_density(aes(x = V1, fill = as.factor(chain_no))) +
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
  inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  data.table::setkey(inf_chain, "i", "sampno", "chain_no")
  n_inf_chain_i <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
  
  if (!is.null(indivs)) {
    use_indivs <- intersect(unique(n_inf_chain_i$i), indivs)
    n_inf_chain_i <- n_inf_chain_i[n_inf_chain_i$i %in% use_indivs, ]
  }
  inf_chain_p_i <- ggplot(n_inf_chain_i) + geom_line(aes(x = sampno, y = V1, col = as.factor(chain_no))) +
    ylab("Estimated total number of infections") +
    xlab("MCMC sample") +
    theme_bw() +
    facet_wrap(~i)
  inf_chain_den_i <- ggplot(n_inf_chain_i) + geom_histogram(aes(x = V1, fill = as.factor(chain_no)), binwidth = 1) +
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
  inf_chain_p <- ggplot(n_inf_chain) + geom_line(aes(x = sampno, y = total_infections, col = as.factor(chain_no))) +
    ylab("Estimated total number of infections") +
    xlab("MCMC sample") +
    theme_bw()
  inf_chain_den <- ggplot(n_inf_chain) + geom_density(aes(x = total_infections, fill = as.factor(chain_no))) +
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
  data.table::setkey(inf_chain, "i", "sampno", "chain_no")
  ## For each individual, how many infections did they have in each sample in total?
  n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]
  
  ## Get quantiles on total number of infections per indiv across all samples
  indiv_hist <- plyr::ddply(n_inf_chain, .(i), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
  colnames(indiv_hist) <- c("individual", "lower", "median", "upper")
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
#' @param burnin only plot samples where sampno > burnin
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
  inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
  inf_chain <- inf_chain[inf_chain$i %in% indivs,]
  if (is.null(inf_chain$chain_no)) {
    inf_chain$chain_no <- 1
  }
  if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)
  
  samps <- sample(unique(inf_chain$sampno), nsamp)
  inf_chain <- inf_chain[inf_chain$sampno %in% samps, ]
  
  ## Get number of probability that infected in a given time point per individual and year
  inf_chain1 <- inf_chain[inf_chain$i %in% indivs, ]
  if (!is.null(subset_times)) inf_chain1 <- inf_chain1[inf_chain1$j %in% subset_times, ]
  data.table::setkey(inf_chain1, "i", "j", "chain_no")
  # max_sampno <- length(unique(inf_chain1$sampno))
  max_sampno <- nsamp
  
  ## Number of samples with a 1 divided by total samples
  densities <- inf_chain1[, list(V1 = sum(x) / max_sampno), by = key(inf_chain1)]
  
  ## Convert to real time points
  densities$j <- as.numeric(possible_exposure_times[densities$j])
  densities$i <- as.numeric(densities$i)
  
  ## If someone wasn't infected in a given year at all, then need a 0
  possible_exposure_times1 <- possible_exposure_times
  if (!is.null(subset_years)) possible_exposure_times1 <- possible_exposure_times[subset_years]
  all_combos <- data.table(expand.grid(i = indivs, j = possible_exposure_times1, chain_no = unique(inf_chain$chain_no)))
  all_combos$j <- as.numeric(all_combos$j)
  all_combos$i <- as.numeric(all_combos$i)
  all_combos <- data.table::fsetdiff(all_combos[, c("i", "j", "chain_no")], densities[, c("i", "j", "chain_no")])
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
  }
  density_plot <- ggplot() +
    geom_line(data = densities, aes(x = j, y = V1, col = as.factor(chain_no)))
  if (!is.null(real_inf_hist)) {
    density_plot <- density_plot +
      geom_vline(data = infection_history1, aes(xintercept = variable), col = "red")
  }
  density_plot <- density_plot +
    facet_wrap(~i) +
    xlab("Time") +
    ylab("Density") +
    theme(legend.position = "bottom") +
    theme_bw()
  
  ## Generate lower, upper and median cumulative infection histories from the
  ## MCMC chain
  tmp_inf_chain <- inf_chain[inf_chain$i %in% indivs, ]
  hist_profiles <- ddply(tmp_inf_chain, .(i, sampno, chain_no), function(x) {
    empty <- numeric(length(possible_exposure_times))
    empty[x[x$x == 1, "j"]] <- 1
    cumsum(empty)
  })
  
  hist_profiles <- hist_profiles[, colnames(hist_profiles) != "sampno"]
  
  colnames(hist_profiles) <- c("i", "chain_no", possible_exposure_times)
  hist_profiles_lower <- ddply(hist_profiles, .(i, chain_no), function(x) apply(x, 2, function(y) quantile(y, c(0.025))))
  hist_profiles_lower <- reshape2::melt(hist_profiles_lower, id.vars = c("i", "chain_no"))
  colnames(hist_profiles_lower) <- c("individual", "chain_no", "variable", "lower")
  
  hist_profiles_upper <- ddply(hist_profiles, .(i, chain_no), function(x) apply(x, 2, function(y) quantile(y, c(0.975))))
  hist_profiles_upper <- reshape2::melt(hist_profiles_upper, id.vars = c("i", "chain_no"))
  colnames(hist_profiles_upper) <- c("individual", "chain_no", "variable", "upper")
  
  hist_profiles_median <- ddply(hist_profiles, .(i, chain_no), function(x) apply(x, 2, function(y) quantile(y, c(0.5))))
  hist_profiles_median <- reshape2::melt(hist_profiles_median, id.vars = c("i", "chain_no"))
  colnames(hist_profiles_median) <- c("individual", "chain_no", "variable", "median")
  
  ## Merge these quantiles into a data frame for plotting
  quant_hist <- merge(hist_profiles_lower, hist_profiles_upper, by = c("individual", "chain_no", "variable"))
  quant_hist <- merge(quant_hist, hist_profiles_median, by = c("individual", "chain_no", "variable"))
  ## If available, process the real infection history matrix for plotting
  real_hist_profiles <- NULL
  if (!is.null(real_inf_hist)) {
    real_hist_profiles <- as.data.frame(t(apply(real_inf_hist, 1, cumsum)))
    colnames(real_hist_profiles) <- possible_exposure_times
    
    real_hist_profiles <- real_hist_profiles[indivs, ]
    real_hist_profiles$individual <- indivs
    real_hist_profiles <- reshape2::melt(real_hist_profiles, id.vars = "individual")
  }
  
  ## Process starting point from MCMC chain
  if (!is.null(start_inf)) {
    start_hist_profiles <- as.data.frame(t(apply(start_inf, 1, cumsum)))
    colnames(start_hist_profiles) <- possible_exposure_times
    start_hist_profiles <- start_hist_profiles[indivs, ]
    start_hist_profiles$individual <- indivs
    start_hist_profiles <- reshape2::melt(start_hist_profiles, id.vars = "individual")
  }
  
  p1 <- ggplot(quant_hist[quant_hist$individual %in% indivs, ]) +
    geom_line(aes(x = as.integer(as.character(variable)), y = median, col = as.factor(chain_no))) +
    geom_ribbon(aes(x = as.integer(as.character(variable)), ymin = lower, ymax = upper, fill = as.factor(chain_no)), alpha = 0.2)
  
  if (!is.null(real_inf_hist)) {
    p1 <- p1 +
      geom_line(data = real_hist_profiles[real_hist_profiles$individual %in% indivs, ], aes(x = as.integer(as.character(variable)), y = value), col = "blue")
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