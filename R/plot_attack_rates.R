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
#' @param plot_ribbon if TRUE, plots a ribbon over time for the attack rate estimates, otherwise plots a pointrange plot
#' @param true_ar data frame of true attack rates, with first column `time` matching `possible_exposure_times`, and second column `AR` giving the attack rate. Column names: population_group, time, AR
#' @param by_group if TRUE, facets the plot by population_group ID
#' @param group_subset if not NULL, plots only this subset of groups eg. 1:5
#' @param plot_residuals if TRUE, plots the residuals between inferred and true attack rate
#' @param colour_by_taken if TRUE, then colours the attack rates by whether or not titres against the circulating antigen at that time were measured
#' @param by_val frequency of x-axis labels
#' @param settings if not NULL, list of serosolver settings as returned from the main serosolver function
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates <- function(infection_histories, 
                              antibody_data=NULL, 
                              demographics = NULL,
                              par_tab=NULL,
                              possible_exposure_times=NULL, 
                              n_alive = NULL,
                              pointsize = 1, fatten = 1,
                              pad_chain = FALSE, prior_pars = NULL,
                              plot_den = FALSE,
                              plot_ribbon=FALSE,
                              true_ar = NULL, by_group = FALSE,
                              group_subset = NULL, plot_residuals = FALSE,
                              colour_by_taken = TRUE, by_val = 5,
                              min_time=min(possible_exposure_times),max_time=max(possible_exposure_times),
                              settings=NULL,
                              verbose=FALSE) {
  ## Some year/sample combinations might have no infections there.
  ## Need to make sure that these get considered
  if (is.null(infection_histories$chain_no)) {
    infection_histories$chain_no <- 1
  }
  
  ## If the list of serosolver settings was included, use these rather than passing each one by one
  if(!is.null(settings)){
    if(verbose) message("Using provided serosolver settings list")
    if(is.null(possible_exposure_times)){
      possible_exposure_times <- settings$possible_exposure_times
      min_time<-min(possible_exposure_times)
      max_time<-max(possible_exposure_times) 
    }
    if(is.null(antibody_data)) antibody_data <- settings$antibody_data
    if(is.null(par_tab)) par_tab <- settings$par_tab
    if(is.null(demographics)) demographics <- settings$demographics
  }
  ## Add stratifying variables to antibody_data and demographics
  ## Setup data vectors and extract
  tmp <- get_demographic_groups(par_tab,antibody_data,demographics, NULL)
  use_demographic_groups <- tmp$use_demographic_groups
  use_timevarying_groups <- tmp$timevarying_demographics
  tmp <- add_stratifying_variables(antibody_data, demographics, par_tab, use_demographic_groups)
  group_ids_vec <- tmp$indiv_pop_group_indices
  
  antibody_data <- tmp$antibody_data
  demographics <- tmp$timevarying_demographics
  demographic_groups <- tmp$demographics
  population_groups <- tmp$population_groups
  if (pad_chain) infection_histories <- pad_inf_chain(infection_histories)
  ## Subset of groups to plot
  if (is.null(group_subset)) {
    group_subset <- unique(antibody_data$population_group)
  }
  
  if (!by_group) {
    antibody_data$population_group <- 1
    infection_histories$population_group <- 1
    if(!is.null(demographics)) demographics$population_group <- 1
    if(!is.null(true_ar)) true_ar$population_group <- 1
  }
  
  
  ## Find inferred total number of infections from the MCMC output
  ## Scale by number of individuals that were alive in each epoch
  ## and generate quantiles
  if (is.null(n_alive)) {
    n_alive <- get_n_alive_group(antibody_data, possible_exposure_times,demographics)
  }
  n_alive <- as.data.frame(n_alive)
  n_alive$population_group <- 1:nrow(n_alive)
  
  unique_groups1 <- unique(antibody_data$population_group)
  n_groups <- length(unique_groups1[!is.na(unique_groups1)])
  n_alive_tot <- get_n_alive(antibody_data, possible_exposure_times)
  infection_histories$j <- possible_exposure_times[infection_histories$j]
  colnames(infection_histories)[1] <- "individual"
  if (!by_group) {
    infection_histories <- merge(infection_histories, data.table(unique(antibody_data[, c("individual", "population_group")])), by = c("individual","population_group"))
  } else {
    if(!is.null(demographics)){
      infection_histories <- merge(infection_histories, 
                                   demographics %>% select(individual, time, population_group) %>% distinct() %>%
                                     rename(j = time) %>%
                                     data.table(), by = c("individual","j"))
    } else { 
      infection_histories <- merge(infection_histories, data.table(unique(antibody_data[, c("individual", "population_group")])), by = c("individual"))
    }
  }
  years <- c(possible_exposure_times, max(possible_exposure_times) + 2)
  data.table::setkey(infection_histories, "samp_no", "j", "chain_no", "population_group")
  tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
  tmp$taken <- tmp$j %in% unique(antibody_data$sample_time)
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
      samp_no = 1:length(prior_dens), j = max(tmp$j) + 1,
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
  n_alive_tmp$j <- possible_exposure_times[n_alive_tmp$j]
  
  tmp <- merge(tmp, data.table(n_alive_tmp), by = c("population_group", "j"))
  tmp <- tmp %>% filter(n_alive > 0)
  tmp$V1 <- tmp$V1 / tmp$n_alive
  
  year_breaks <- c(min_time, seq(by_val * round(min_time / by_val), max_time, by = by_val))
  year_labels <- c(min_time, seq(by_val * round(min_time / by_val), max_time, by = by_val))
  
  if (!is.null(prior_dens)) {
    year_breaks <- c(year_breaks, max_time + 2)
    year_labels <- c(year_labels, "Prior")
  }
  if (!plot_den) {
    quantiles <- tmp %>% group_by(population_group,j) %>% dplyr::summarise(lower=quantile(V1, 0.025),median=quantile(V1,0.5),upper=quantile(V1,0.975)) %>% ungroup()
    quantiles$taken <- quantiles$j %in% unique(antibody_data$sample_time)
    quantiles$taken <- ifelse(quantiles$taken, "Yes", "No")
    
    quantiles$tested <- quantiles$j %in% unique(antibody_data$biomarker_id)
    quantiles$tested <- ifelse(quantiles$tested, "Yes", "No")
    
    if (!is.null(prior_dens)) {
      quantiles[quantiles$j == max(years), "taken"] <- "Prior"
    }
    
    colnames(quantiles)[which(colnames(quantiles) == "taken")] <- "Sample taken"
    colnames(quantiles)[which(colnames(quantiles) == "tested")]  <- "Biomarker tested"
    quantiles$population_group <- as.factor(quantiles$population_group)
    p <- ggplot(quantiles[quantiles$population_group %in% group_subset, ]) +
      coord_cartesian(ylim=c(0,1)) +
      theme_classic() +
      ylab("Estimated attack rate") +
      xlab("Year")
    
    ## Colour depending on whether or not titres were taken in each year
    if(!plot_ribbon){
      if (colour_by_taken == TRUE) {
        p <- p + geom_pointrange(aes(
          x = j, y = median, ymin = lower, ymax = upper,
          col = `Sample taken`
        ),
        size = pointsize,
        fatten = fatten
        )  +
          scale_fill_manual(name="Samples taken",values=c("No"="darkorange","Yes"="blue","Prior"="grey40"))+
          scale_color_manual(name="Samples taken",values=c("No"="darkorange","Yes"="blue","Prior"="grey40"))
      } else {
        p <- p + geom_pointrange(aes(
          x = j, y = median, ymin = lower, ymax = upper,
          col = `Biomarker tested`
        ),
        size = pointsize,
        fatten = fatten
        )+
          scale_fill_manual(name="Biomarker tested",values=c("No"="darkorange","Yes"="blue","Prior"="grey40"))+
          scale_color_manual(name="Biomarker tested",values=c("No"="darkorange","Yes"="blue","Prior"="grey40"))
      }
    } else {
      if(!by_group){
        p <- p + geom_ribbon(aes(x = j, ymin = lower, ymax = upper), alpha = 0.25) +
          geom_line(aes(x = j, y = median), linewidth = 0.75)
      } else {
        p <- p + geom_ribbon(aes(x = j, ymin = lower, ymax = upper,fill=population_group,group=population_group), alpha = 0.25) +  geom_line(aes(x = j, y = median,col=population_group,group=population_group), linewidth = 0.75) +
          scale_fill_viridis_d(name="Population group") +
          scale_color_viridis_d(name="Population group")
      }
    }
  } else {
    #tmp$j <- years[tmp$j]
    colnames(tmp)[which(colnames(tmp) == "j")] <- "time"
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
      geom_violin(aes(x = time, y = V1, fill = taken, group = time),
                  alpha=0.25,
                  draw_quantiles = c(0.025,0.5,0.975), scale = "width",
                  adjust=2
      )+
      scale_fill_manual(name="Samples taken",values=c("No"="darkorange","Yes"="blue","Prior"="grey40"))+
      scale_color_manual(name="Samples taken",values=c("No"="darkorange","Yes"="blue","Prior"="grey40"))
  }
  if (!is.null(true_ar) & !plot_residuals) {
    p <- p +
      geom_point(
        data = true_ar[true_ar$population_group %in% group_subset, ], aes(x = time, y = AR,shape="True attack rate"),stroke=1.25,
        col = "black", size = 2.5
      ) +
      scale_shape_manual(name="",values=c("True attack rate"=1))
  }
  if (!plot_residuals) {
    p <- p +
      scale_y_continuous(expand = c(0, 0)) 
  }
  
  if (by_group & !plot_ribbon) {
    p <- p + facet_wrap(~population_group, ncol = 1)
  }
  if (!is.null(prior_dens)) {
    max_time <- max_time + 2.5
  }
  
  p <- p +
    scale_x_continuous(breaks = year_breaks, labels = year_labels,limits=c(min_time-0.5,max_time+0.5)) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ylab("Estimated attack rate") +
    xlab("Date")
  
  if (plot_residuals) {
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))
    p_res <- ggplot(tmp) + geom_density(aes(x = V1), fill = "grey40",alpha=0.5) +
      facet_wrap(~population_group) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "blue") +
      scale_x_continuous(limits = c(-1, 1)) +
      theme_classic() + xlab("Estimated AR - true AR") + ylab("Density")
    return(list(p, p_res))
  }
  
  return(p)
}