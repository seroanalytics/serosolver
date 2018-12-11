#' Formatted quantiles
#'
#' Given a vector of MCMC samples, generates and formats the desired quantile estimates
#' @param x the vector to summarise
#' @param sig_f how many significant figures to print
#' @param qs the vector of quantiles
#' @param as_text if TRUE, formats nicely as text rather than a vector of numbers
#' @return the formatted quantiles
#' @export
generate_quantiles <- function(x, sig_f = 3, qs = c(0.025, 0.5, 0.975), as_text = TRUE) {
  res <- signif(quantile(x, qs), sig_f)
  if (as_text) {
    res <- paste(res[2], " (", res[1], "-", res[3], ")", sep = "")
  }
  return(res)
}


#' Generate titre credible intervals
#'
#' Generates credible intervals on titres and infection histories from an MCMC chain output.
#' @param chain the full MCMC chain to generate titre trajectories from
#' @param infection_histories the MCMC chain for infection histories
#' @param titre_dat the data frame of titre data
#' @param individuals the subset of individuals to generate credible intervals for
#' @param antigenic_map the unmelted antigenic map
#' @param par_tab the table controlling the parameters in the MCMC chain
#' @param nsamp number of samples to take from posterior
#' @param add_residuals if true, returns an extra output summarising residuals between the model prediction and data
#' @param mu_indices vector of integers. for random effects on boosting parameter, mu. If random mus are included in the parameter table, this vector specifies which mu to use for each circulation year. For example, if years 1970-1976 have unique boosting, then mu_indices should be c(1,2,3,4,5,6). If every 3 year block shares has a unique boosting parameter, then this should be c(1,1,1,2,2,2)
#' @param measurement_indices_by_time default NULL, optional vector giving the index of `measurement_bias` that each strain uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clusters
#' @param for_res_plot TRUE/FALSE value. If using the output of this for plotting of residuals, returns the actual data points rather than summary statistics
#' @return a list with the titre predictions (95% credible intervals, median and multivariate posterior mode) and the probabilities of infection for each individual in each epoch
#' @export
get_titre_predictions <- function(chain, infection_histories, titre_dat,
                                  individuals, antigenic_map,
                                  par_tab,
                                  nsamp = 100, add_residuals = FALSE,
                                  mu_indices = NULL,
                                  measurement_indices_by_time = NULL,
                                  for_res_plot = FALSE) {
  ## Need to align the iterations of the two MCMC chains
  ## and choose some random samples
  samps <- intersect(unique(infection_histories$sampno), unique(chain$sampno))
  chain <- chain[chain$sampno %in% samps, ]
  infection_histories <- infection_histories[infection_histories$sampno %in% samps, ]

  ## Take subset of individuals
  titre_dat <- titre_dat[titre_dat$individual %in% individuals, ]
  infection_histories <- infection_histories[infection_histories$i %in% individuals, ]

  titre_dat$individual <- match(titre_dat$individual, individuals)
  infection_histories$i <- match(infection_histories$i, individuals)

  ## Format the antigenic map to solve the model
  strain_isolation_times <- unique(antigenic_map$inf_years)
  nstrain <- length(strain_isolation_times)
  n_indiv <- length(individuals)

  ## Empty data structures to save output to
  infection_history_dens <- NULL
  tmp_samp <- sample(samps, nsamp)

  ## See the function in posteriors.R
  f <- create_posterior_func(par_tab, titre_dat, antigenic_map, 100, mu_indices = mu_indices,
                             measurement_indices_by_time = measurement_indices_by_time, function_type = 3)

  predicted_titres <- residuals <- matrix(nrow = nrow(titre_dat), ncol = nsamp)
  samp_record <- numeric(nsamp)

  ## For each sample, take values for theta and infection histories and simulate titres
  for (i in 1:nsamp) {
    index <- tmp_samp[i]
    pars <- get_index_pars(chain, which(chain$sampno == index))
    tmp_inf_hist <- infection_histories[infection_histories$sampno == index, ]
    tmp_inf_hist <- as.matrix(Matrix::sparseMatrix(i = tmp_inf_hist$i, j = tmp_inf_hist$j, x = tmp_inf_hist$x, dims = c(n_indiv, nstrain)))
    predicted_titres[, i] <- f(pars, tmp_inf_hist)

    ## Get residuals between observations and predictions
    residuals[, i] <- titre_dat$titre - floor(predicted_titres[, i])
    samp_record[i] <- index
  }

  colnames(predicted_titres) <- tmp_samp

  ## If generating for residual plot, can return now
  if (for_res_plot) return(list(residuals, samp_record, titre_dat, predicted_titres))

  residuals <- cbind(titre_dat, residuals)
  ## Get 95% credible interval and means
  dat2 <- t(apply(predicted_titres, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
  residuals <- t(apply(residuals, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
  residuals <- cbind(titre_dat, residuals)

  ## Find multivariate posterior mode estimate from the chain
  best_pars <- get_best_pars(chain)
  best_I <- chain$sampno[which.max(chain$lnlike)]
  best_inf <- infection_histories[infection_histories$sampno == best_I, ]
  best_inf <- as.matrix(Matrix::sparseMatrix(i = best_inf$i, j = best_inf$j, x = best_inf$x, dims = c(n_indiv, nstrain)))

  ## Generate trajectory for best parameters
  best_traj <- f(best_pars, best_inf)
  best_residuals <- titre_dat$titre - floor(best_traj)
  best_residuals <- cbind(titre_dat, best_residuals, "sampno" = best_I)
  dat2 <- as.data.frame(dat2)
  colnames(dat2) <- c("lower", "median", "upper")
  dat2$max <- best_traj
  dat2[dat2 < 0] <- 0
  dat2 <- cbind(titre_dat, dat2)
  tmp_inf_chain <- data.table(subset(infection_histories, sampno %in% tmp_samp))

  ## Get infection history density for each individual and each epoch
  data.table::setkey(tmp_inf_chain, "i", "j")
  infection_history_dens <- tmp_inf_chain[, list(V1 = sum(x) / length(tmp_samp)), by = key(tmp_inf_chain)]
  infection_history_dens$j <- strain_isolation_times[infection_history_dens$j]
  colnames(infection_history_dens) <- c("individual", "variable", "value")
  infection_history_final <- NULL

  ## For each individual, get density for the probability that an epoch was an infection time
  ## The point of the following loop is to mask the densities where infection epochs were either
  ## before an individual was born or after the time that a blood sample was taken
  for (indiv in unique(infection_history_dens$individual)) {
    sample_times <- unique(titre_dat[titre_dat$individual == indiv, "samples"])
    tmp <- NULL
    for (samp in sample_times) {
      indiv_inf_hist <- infection_history_dens[infection_history_dens$individual == indiv, ]
      indiv_inf_hist[indiv_inf_hist$variable > samp, "value"] <- 0
      indiv_inf_hist <- cbind(indiv_inf_hist, "samples" = samp)
      tmp <- rbind(tmp, indiv_inf_hist)
    }
    infection_history_final <- rbind(infection_history_final, tmp)
  }
  dat2$individual <- individuals[dat2$individual]
  infection_history_final$individual <- individuals[infection_history_final$individual]
  if (add_residuals) {
    result <- list("predictions" = dat2, "histories" = infection_history_final, "residuals" = residuals, "bestRes" = best_residuals)
  } else {
    result <- list("predictions" = dat2, "histories" = infection_history_final)
  }
  return(result)
}

#' Plots infection histories and titre model fits
#'
#' Given outputs from an MCMC run and the data used for fitting, generates an NxM matrix of plots where N is the number of individuals to be plotted and M is the range of sampling times. Where data are available, plots the observed titres and model predicted trajectories
#' @inheritParams get_titre_predictions
#' @return a ggplot2 object
#' @export
plot_infection_histories <- function(chain, infection_histories, titre_dat,
                                     individuals, antigenic_map, par_tab,
                                     nsamp = 100,
                                     mu_indices = NULL,
                                     measurement_indices_by_time = NULL) {
  individuals <- individuals[order(individuals)]

  ## Generate titre predictions
  tmp <- get_titre_predictions(
    chain, infection_histories, titre_dat, individuals,
    antigenic_map, par_tab, nsamp, FALSE, mu_indices,
    measurement_indices_by_time
  )

  ## Use these titre predictions and summary statistics on infection histories
  dens <- tmp[[1]]
  infection_history <- tmp[[2]]
  p <- ggplot(dens) +
    geom_line(aes(x = virus, y = median), col = "blue") +
    geom_ribbon(aes(x = virus, ymin = lower, ymax = upper), alpha = 0.25, fill = "blue") +
    geom_vline(data = infection_history, aes(xintercept = variable, alpha = value)) +
    geom_point(data = dens, aes(x = virus, y = titre), col = "red", size = 0.5) +
    facet_grid(individual ~ samples) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    scale_alpha(limits = c(0, 1), range = c(0, 1)) +
    xlab("Year") +
    ylab("Titre") +
    scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))
  p
}

#' Plot inferred posterior distributions
#'
#' Produces and saves estimated posterior distributions
#' @param chain the full MCMC chain to generate titre trajectories from
#' @param par_tab the table controlling the parameters in the MCMC chain
#' @param burnin if not already discarded, discard burn in from chain (takes rows where sampno > burnin)
#' @param samples how many samples from the chain to take
#' @param calculate_ess if TRUE, calculates the ESS for all free parameters
#' @param plot_corr if TRUE, returns a pairwise correlation plot of free parameters
#' @param save_plots if TRUE, directly saves the plots as svgs
#' @param plot_mcmc if TRUE, plots the MCMC chain traces
#' @param save_loc the full directory path of where to save plots
#' @return a list of all ggplot objects
#' @export
plot_posteriors <- function(chain,
                            par_tab,
                            burnin = 0,
                            samples = 1000,
                            calculate_ess = TRUE,
                            plot_corr = TRUE,
                            save_plots = FALSE,
                            plot_mcmc = TRUE,
                            save_loc = "") {
  ## Combined chain
  free_chain <- chain[chain$sampno > burnin, ]
  free_chain <- chain[, par_tab[par_tab$fixed == 0, "names"]]
  thin_free_chain <- free_chain[sample(1:nrow(free_chain), samples, replace = T), ]

  ## Plot correlations
  corP <- NULL
  if (plot_corr) {
    corP <- GGally::ggpairs(thin_free_chain) + theme_bw()

    if (save_plots) {
      to.svg(print(corP), paste0(save_loc, "corPlot.svg"))
      to.svg(coda::autocorr.plot(free_chain), paste0(save_loc, "autocorr.svg"))
      # to.pdf(print(corP), paste0(save_loc,"corPlot.pdf"))
      # to.pdf(coda::autocorr.plot(free_chain),paste0(save_loc,"autocorr.pdf"))
    }
  }

  ## Get quantiles and create table of results
  ## These are from Adam's summaries...
  thin_free_chain$sigma1drop <- thin_free_chain$mu * thin_free_chain$sigma1
  thin_free_chain$sigma2drop <- thin_free_chain$mu_short * thin_free_chain$sigma2
  thin_free_chain$wane_yr <- thin_free_chain$mu_short * thin_free_chain$wane
  thin_free_chain$errorCorrect1 <- pnorm(3, mean = 1.5, sd = thin_free_chain$error) - pnorm(2, mean = 2.5, sd = thin_free_chain$error)
  thin_free_chain$errorCorrect2 <- pnorm(4, mean = 1.5, sd = thin_free_chain$error) - pnorm(1, mean = 2.5, sd = thin_free_chain$error)

  ## Calculate ESS for all free parameters
  all_ess <- NULL
  if (calculate_ess) {
    all_ess <- coda::effectiveSize(free_chain)
    all_ess <- c(all_ess, rep(NA, 5))
  }

  densP <- traceP <- NULL
  if (plot_mcmc) {
    densP <- bayesplot::mcmc_hist(free_chain) + theme()
    traceP <- bayesplot::mcmc_trace(free_chain) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    if (save_plots) {
      to.svg(print(densP), paste0(save_loc, "densities.svg"))
      to.svg(print(traceP), paste0(save_loc, "traces.svg"))
      # to.pdf(print(densP), paste0(save_loc,"densities.pdf"))
      # to.pdf(print(traceP), paste0(save_loc,"traces.pdf"))
    }
  }

  results <- apply(thin_free_chain, 2, function(x) generate_quantiles(x))
  all_results <- data.frame(
    parameter = c(
      par_tab[which(par_tab$fixed == 0), "names"],
      "sigma1drop", "sigma2drop", "wane_yr", "errorCorrect1", "errorCorrect2"
    ),
    estimate = results,
    ESS = all_ess
  )
  if (save_plots) write.table(all_results, paste0(save_loc, "estimates.csv"), row.names = FALSE, sep = ",")
  return(list(results = all_results, corPlot = corP, densP = densP, traceP = traceP))
}

#' Get posterior information infection histories
#'
#' Finds the median, mean and 95% credible intervals for the attack rates and total number of infections per individual
#' @inheritParams plot_infection_history_chains_time
#' @param years vector of years giving the labels to matching j in the infection_histories data frame
#' @return a list of data frames with summary statistics
#' @export
calculate_infection_history_statistics <- function(inf_chain, burnin=0, years=NULL, n_alive=NULL){
    inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
    inf_chain <- pad_inf_chain(inf_chain)

      
    data.table::setkey(inf_chain, "j", "sampno")
    n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]

    if (!is.null(n_alive)) {
        n_inf_chain$V1 <- n_inf_chain$V1 / n_alive[n_inf_chain$j]
    }
     if (!is.null(years)) {
         n_inf_chain$j <- years[n_inf_chain$j]
     }

    data.table::setkey(n_inf_chain, "j")
    n_inf_chain_summaries <- n_inf_chain[,list(mean=mean(V1),median=median(V1),
                                               lower_quantile=quantile(V1,c(0.025)),
                                               upper_quantile=quantile(V1,c(0.975))),
                                         by=key(n_inf_chain)]

    

    data.table::setkey(inf_chain, "i", "sampno")
    n_inf_chain_i <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]

    data.table::setkey(n_inf_chain_i, "i")
    n_inf_chain_i_summaries <- n_inf_chain_i[,list(mean=mean(V1),median=median(V1),
                                               lower_quantile=quantile(V1,c(0.025)),
                                               upper_quantile=quantile(V1,c(0.975))),
                                         by=key(n_inf_chain_i)]

    return(list("by_year"=n_inf_chain_summaries,"by_indiv"=n_inf_chain_i_summaries))    
}

#' Plot historical attack rates monthly
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories for monthly. The main difference compared to the normal attack rate plot is that pointrange plots don't make as much sense at a very fine time resolution.
#' @inheritParams plot_attack_rates
#' @param ymax Numeric. the maximum y value to put on the axis
#' @param buckets Integer. How many buckets of time is each year split into? ie. 12 for monthly data, 4 for quarterly etc.
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates_monthly <- function(infection_histories, titre_dat, years, n_alive = NULL, ymax = 0.1, buckets = 1) {
  ## Some year/sample combinations might have no infections there.
  ## Need to make sure that these get considered
  infection_histories <- pad_inf_chain(infection_histories)
  ## Find inferred total number of infections from the MCMC output
  ## Scale by number of individuals that were alive in each epoch
  ## and generate quantiles
  months <- 1:buckets
  years <- range(floor(years / buckets))
  years <- years[1]:years[2]
  labels <- c(sapply(years, function(x) paste0(months, "/", x)))
  labels1 <- labels[1:length(years)]
  labels1 <- labels1[seq(1, length(labels1), by = buckets)]
  year_break <- years[seq(1, length(years), by = buckets)]
  if (is.null(n_alive)) {
    n_alive <- get_n_alive(titre_dat, years)
  }

  ## Sum infections per year for each MCMC sample
  data.table::setkey(infection_histories, "sampno", "j")
  tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]

  quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
  colnames(quantiles) <- c("j", "lower", "median", "upper")
  quantiles[c("lower", "median", "upper")] <- quantiles[c("lower", "median", "upper")] / n_alive
  quantiles$year <- years[quantiles$j]
  quantiles$taken <- quantiles$year %in% unique(titre_dat$samples)

  ## Colour depending on whether or not titres were taken in each year
  quantiles$taken <- ifelse(quantiles$taken, "Yes", "No")

  p <- ggplot(quantiles) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
    geom_line(aes(x = year, y = median), col = "red") +
    geom_point(aes(x = year, y = median), col = "purple", size = 0.5) +
    scale_y_continuous(limits = c(-0.005, ymax), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), breaks = year_break, labels = labels1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Estimated monthly per capita incidence") +
    xlab("Date")
  return(p)
}

#' Plot historical attack rates
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories
#' @param infection_histories the MCMC chain for infection histories
#' @param titre_dat the data frame of titre data
#' @param year vector of the epochs of potential circulation
#' @param n_alive vector with the number of people alive in each year of circulation. Can be left as NULL, and ages will be used to infer this
#' @param pointsize Numeric - how big should each point be?
#' @param fatten Numeric - fatten parameter for ggplot pointrange
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates <- function(infection_histories, titre_dat, years, n_alive = NULL, pointsize = 1, fatten = 1) {
  ## Some year/sample combinations might have no infections there.
  ## Need to make sure that these get considered
  infection_histories <- pad_inf_chain(infection_histories)
  ## Find inferred total number of infections from the MCMC output
  ## Scale by number of individuals that were alive in each epoch
  ## and generate quantiles
  if (is.null(n_alive)) {
    n_alive <- get_n_alive(titre_dat, years)
  }

  data.table::setkey(infection_histories, "sampno", "j")
  tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
  quantiles <- ddply(tmp, ~j, function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
  colnames(quantiles) <- c("j", "lower", "median", "upper")
  quantiles[c("lower", "median", "upper")] <- quantiles[c("lower", "median", "upper")] / n_alive
  quantiles$year <- years[quantiles$j]
  quantiles$taken <- quantiles$year %in% unique(titre_dat$samples)

  ## Colour depending on whether or not titres were taken in each year
  quantiles$taken <- ifelse(quantiles$taken, "Yes", "No")

  p <- ggplot(quantiles) +
    geom_pointrange(aes(x = year, y = median, ymin = lower, ymax = upper, col = taken), size = pointsize, fatten = fatten) +
    scale_y_continuous(limits = c(-0.1, 1), expand = c(0, 0)) +
    theme_bw() +
    ylab("Estimated attack rate") +
    xlab("Year")
  return(p)
}

#' Pad infection history chain
#'
#' Given that the MCMC sampler only stores present infections (ie. there are no entries for 0s from the infection history matrix), for some summaries we need to add these 0s back in to avoid bias.
#' @param inf_chain the data table with infection history samples from \code{\link{run_MCMC}}
#' @return the same inf_chain that was passed in, but with 0s for missing i/j/sampno combinations
#' @export
pad_inf_chain <- function(inf_chain) {
  expanded_values <- data.table(expand.grid(
    i = unique(inf_chain$i),
    j = unique(inf_chain$j),
    sampno = unique(inf_chain$sampno)
  ))
  diff_infs <- fsetdiff(expanded_values, inf_chain[, c("i", "j", "sampno")])
  diff_infs$x <- 0
  inf_chain <- rbind(inf_chain, diff_infs)
  return(inf_chain)
}

#' Plot MCMC trace for infections per year
#'
#' @param inf_chain the data table with infection history samples from \code{\link{run_MCMC}}
#' @param burnin optionally remove all sampno < burnin from the chain
#' @param years vector of integers, if not NULL, only plots a subset of years (where 1 is the first year eg. 1968)
#' @param n_alive if not NULL, then divides number of infections per year by number alive to give attack rates rather than total infections
#' @return a list of two ggplot objects - the MCMC trace and MCMC densities
#' @seealso \code{\link{plot_infection_history_chains_indiv}}
#' @export
plot_infection_history_chains_time <- function(inf_chain, burnin = 0, years = NULL, n_alive = NULL) {
  inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
  inf_chain <- pad_inf_chain(inf_chain)
  data.table::setkey(inf_chain, "j", "sampno")
  n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]

  if (!is.null(n_alive)) {
    n_inf_chain$V1 <- n_inf_chain$V1 / n_alive[n_inf_chain$j]
  }

  if (!is.null(years)) {
    use_years <- intersect(unique(n_inf_chain$j), years)
    n_inf_chain <- n_inf_chain[n_inf_chain$j %in% use_years, ]
  }

  inf_chain_p <- ggplot(n_inf_chain) + geom_line(aes(x = sampno, y = V1)) + facet_wrap(~j)
  inf_chain_den <- ggplot(n_inf_chain) + geom_density(aes(x = V1)) + facet_wrap(~j)
  return(list(inf_chain_p, inf_chain_den))
}

#' Plot MCMC trace for infections per individual
#'
#' @inheritParams plot_infection_history_chains_time
#' @param indivs vector of integers, if not NULL, only plots a subset of individuals (where 1 is the first individual)
#' @return a list of two ggplot objects - the MCMC trace and MCMC densities
#' @seealso \code{\link{plot_infection_history_chains_indiv}}
#' @export
plot_infection_history_chains_indiv <- function(inf_chain, burnin = 0, indivs = NULL) {
  inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
  inf_chain <- pad_inf_chain(inf_chain)
  data.table::setkey(inf_chain, "i", "sampno")
  n_inf_chain_i <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]

  if (!is.null(indivs)) {
    use_indivs <- intersect(unique(n_inf_chain_i$i), indivs)
    n_inf_chain_i <- n_inf_chain_i[n_inf_chain_i$i %in% use_indivs, ]
  }

  inf_chain_p_i <- ggplot(n_inf_chain_i) + geom_line(aes(x = sampno, y = V1)) + facet_wrap(~i)
  inf_chain_den_i <- ggplot(n_inf_chain_i) + geom_density(aes(x = V1)) + facet_wrap(~j)
  return(list(inf_chain_p_i, inf_chain_den_i))
}

#' Plot point range number infections per individual
#'
#' @param inf_chain the data table with infection history samples from \code{\link{run_MCMC}}
#' @return a ggplot object
#' @export
plot_number_infections <- function(inf_chain) {
  n_strain <- max(inf_chain$j)
  data.table::setkey(inf_chain, "i", "sampno")
  ## For each individual, how many infections did they have in each sample in total?
  n_inf_chain <- inf_chain[, list(V1 = sum(x)), by = key(inf_chain)]

  ## Get quantiles on total number of infections per indiv across all samples
  indiv_hist <- plyr::ddply(n_inf_chain, .(i), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
  colnames(indiv_hist) <- c("individual", "lower", "median", "upper")
  indiv_hist <- indiv_hist[order(indiv_hist$median), ]
  indiv_hist$individual <- 1:nrow(indiv_hist)
  p <- ggplot(indiv_hist) +
    geom_pointrange(aes(x = individual, y = median, ymin = lower, ymax = upper),
      size = 0.1, shape = 21, fatten = 0.1
    ) +
    scale_y_continuous(limits = c(0, n_strain), expand = c(0, 0)) +
    xlab("Individual") +
    ylab("Estimated number of infections") +
    theme_bw()
  return(p)
}

#' Useful plot for looking at simulated data
#'
#' Plots measured titres and known infection histories for all individuals, facetted by sample time
#' @param titre_dat the data frame of titre data
#' @param infection_histories the infection history matrix
#' @param strain_isolation_times the vector of times at which individuals could be infected
#' @param n_samps how many individuals to plot
#' @param start_inf if not NULL, plots the infection history matrix used as the starting point in the MCMC chain
#' @return a ggplot object
#' @export
plot_data <- function(titre_dat, infection_histories, strain_isolation_times, n_samps, start_inf = NULL) {
  indivs <- unique(titre_dat$individual)
  infection_history <- as.data.frame(cbind(indivs, infection_histories))
  colnames(infection_history) <- c("individual", strain_isolation_times)
  melted_inf_hist <- reshape2::melt(infection_history, id.vars = "individual")
  melted_inf_hist$variable <- as.numeric(as.character(melted_inf_hist$variable))
  melted_inf_hist <- melted_inf_hist[melted_inf_hist$value > 0, ]
  samps <- sample(unique(titre_dat$individual), n_samps)

  p1 <- ggplot(titre_dat[titre_dat$individual %in% samps, ]) +
    geom_line(aes(x = as.integer(virus), y = titre)) +
    geom_vline(data = melted_inf_hist[melted_inf_hist$individual %in% samps, ], aes(xintercept = variable), col = "red", linetype = "dashed") +
    theme_bw()

  if (!is.null(start_inf)) {
    start_inf_hist <- as.data.frame(cbind(indivs, start_inf))
    colnames(start_inf_hist) <- c("individual", strain_isolation_times)
    melted_start_hist <- reshape2::melt(start_inf_hist, id.vars = "individual")
    melted_start_hist$variable <- as.numeric(as.character(melted_start_hist$variable))
    melted_start_hist <- melted_start_hist[melted_start_hist$value > 0, ]
    p1 <- p1 + geom_vline(data = melted_start_hist[melted_start_hist$individual %in% samps, ], aes(xintercept = variable), col = "blue", linetype = "dashed")
  }
  p1 <- p1 +
    facet_grid(individual ~ samples)
  return(p1)
}

#' Plot cumulative and per time posterior infection probability densities
#'
#' For each individual requested, plots the median and 95% quantiles on a) the cumulative number of infections over a lifetime and b) the posterior probability that an infection occured in a given time point
#' @param inf_chain_file the relative file path to the infection history chain
#' @param burnin only plot samples where sampno > burnin
#' @param indivs vector of individual ids to plot
#' @param real_inf_hist if not NULL, adds lines to the plots showing the known true infection times
#' @param start_inf if not NULL, adds lines to show where the MCMC chain started
#' @param strain_isolation_times vector of times at which individuals could have been infected
#' @param nsamp how many samples from the MCMC chain to take?
#' @param ages if not NULL, adds lines to show when an individual was born
#' @param number_col how many columns to use for the cumulative infection history plot
#' @return two ggplot objects
#' @export
generate_cumulative_inf_plots <- function(inf_chain_file, burnin = 0, indivs, real_inf_hist = NULL, start_inf = NULL,
                                          strain_isolation_times, nsamp = 100, ages = NULL, number_col = 1) {
  inf_chain <- data.table::fread(inf_chain_file)
  inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
  inf_chain <- pad_inf_chain(inf_chain)

  samps <- sample(unique(inf_chain$sampno), nsamp)
  inf_chain <- inf_chain[inf_chain$sampno %in% samps, ]

  ## Get number of probability that infected in a given time point per individual and year
  inf_chain1 <- inf_chain[inf_chain$i %in% indivs, ]
  data.table::setkey(inf_chain1, "i", "j")
  max_sampno <- length(unique(inf_chain1$sampno))

  ## Number of samples with a 1 divided by total samples
  densities <- inf_chain1[, list(V1 = sum(x) / max_sampno), by = key(inf_chain1)]

  ## Convert to real time points
  densities$j <- as.numeric(strain_isolation_times[densities$j])
  densities$i <- as.numeric(densities$i)

  ## If someone wasn't infected in a given year at all, then need a 0
  all_combos <- data.table(expand.grid(i = indivs, j = strain_isolation_times))
  all_combos$j <- as.numeric(all_combos$j)
  all_combos$i <- as.numeric(all_combos$i)
  all_combos <- data.table::fsetdiff(all_combos[, c("i", "j")], densities[, c("i", "j")])
  all_combos$V1 <- 0
  densities <- rbind(all_combos, densities)

  if (!is.null(real_inf_hist)) {
    infection_history1 <- as.data.frame(real_inf_hist)
    infection_history1 <- infection_history1[indivs, ]
    infection_history1$individual <- indivs

    colnames(infection_history1) <- c(strain_isolation_times, "i")
    infection_history1 <- reshape2::melt(infection_history1, id.vars = "i")
    infection_history1$variable <- as.numeric(as.character(infection_history1$variable))
    infection_history1 <- infection_history1[infection_history1$value == 1, ]
  }
  density_plot <- ggplot() +
    geom_line(data = densities, aes(x = j, y = V1))
  if (!is.null(real_inf_hist)) {
    density_plot <- density_plot +
      geom_vline(data = infection_history1, aes(xintercept = variable), col = "red")
  }
  density_plot <- density_plot +
    facet_wrap(~i) +
    xlab("Year") +
    ylab("Density") +
    theme_bw()

  ## Generate lower, upper and median cumulative infection histories from the
  ## MCMC chain
  tmp_inf_chain <- inf_chain[inf_chain$i %in% indivs, ]
  hist_profiles <- ddply(tmp_inf_chain, .(i, sampno), function(x) {
    empty <- numeric(length(strain_isolation_times))
    empty[x$j] <- 1
    cumsum(empty)
  })


  hist_profiles <- hist_profiles[, colnames(hist_profiles) != "sampno"]
  colnames(hist_profiles) <- c("i", strain_isolation_times)
  hist_profiles_lower <- ddply(hist_profiles, ~i, function(x) apply(x, 2, function(y) quantile(y, c(0.025))))
  hist_profiles_lower <- reshape2::melt(hist_profiles_lower, id.vars = "i")
  colnames(hist_profiles_lower) <- c("individual", "variable", "lower")

  hist_profiles_upper <- ddply(hist_profiles, ~i, function(x) apply(x, 2, function(y) quantile(y, c(0.975))))
  hist_profiles_upper <- reshape2::melt(hist_profiles_upper, id.vars = "i")
  colnames(hist_profiles_upper) <- c("individual", "variable", "upper")

  hist_profiles_median <- ddply(hist_profiles, ~i, function(x) apply(x, 2, function(y) quantile(y, c(0.5))))
  hist_profiles_median <- reshape2::melt(hist_profiles_median, id.vars = "i")
  colnames(hist_profiles_median) <- c("individual", "variable", "median")

  ## Merge these quantiles into a data frame for plotting
  quant_hist <- merge(hist_profiles_lower, hist_profiles_upper, by = c("individual", "variable"))
  quant_hist <- merge(quant_hist, hist_profiles_median, by = c("individual", "variable"))

  ## If available, process the real infection history matrix for plotting
  if (!is.null(real_inf_hist)) {
    real_hist_profiles <- as.data.frame(t(apply(real_inf_hist, 1, cumsum)))
    colnames(real_hist_profiles) <- strain_isolation_times

    real_hist_profiles <- real_hist_profiles[indivs, ]
    real_hist_profiles$individual <- indivs
    real_hist_profiles <- reshape2::melt(real_hist_profiles, id.vars = "individual")
  }

  ## Process starting point from MCMC chain
  if (!is.null(start_inf)) {
    start_hist_profiles <- as.data.frame(t(apply(start_inf, 1, cumsum)))
    colnames(start_hist_profiles) <- strain_isolation_times
    start_hist_profiles <- start_hist_profiles[indivs, ]
    start_hist_profiles$individual <- indivs
    start_hist_profiles <- reshape2::melt(start_hist_profiles, id.vars = "individual")
  }

  p1 <- ggplot(quant_hist[quant_hist$individual %in% indivs, ]) +
    geom_line(aes(x = as.integer(as.character(variable)), y = median)) +
    geom_ribbon(aes(x = as.integer(as.character(variable)), ymin = lower, ymax = upper), alpha = 0.2)

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
    age_mask <- create_age_mask(tmp_age[, 2], strain_isolation_times)
    age_dat <- data.frame(j = age_mask, individual = indivs[order(indivs)])
    p1 <- p1 + geom_vline(data = age_dat, aes(xintercept = strain_isolation_times[j]), col = "purple", linetype = "dashed")
  }
  p1 <- p1 + facet_wrap(~individual, ncol = number_col) +
    theme_bw() +
    ylab("Cumulative infections") +
    xlab("Circulation time")
  return(list(p1, density_plot))
}

#' Titre dependent boosting relationship
#'
#' Calculates the inferred titre dependent boosting relationship from the MCMC chain
#' @param chain the MCMC chain
#' @param n number of samples to take
#' @param titres the vector of titres to calculate boosting values at
#' @return a data frame of quantiles for the inferred boost from different titre levels
#' @export
titre_dependent_boosting_plot <- function(chain, n, titres = seq(0, 8, by = 0.1)) {
  sampnos <- sample(unique(chain$sampno), n)
  store <- matrix(nrow = n, ncol = length(titres))
  i <- 1
  for (samp in sampnos) {
    pars <- as.numeric(chain[chain$sampno == samp, ])
    names(pars) <- colnames(chain)
    mu <- pars["mu"] + pars["mu_short"]

    gradient <- pars["gradient"]
    boost_limit <- pars["boost_limit"]
    boost <- mu * (1 - gradient * titres)
    boost[which(titres > boost_limit)] <- mu * (1 - gradient * boost_limit)
    store[i, ] <- boost
    i <- i + 1
  }
  range <- apply(store, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
  return(range)
}
