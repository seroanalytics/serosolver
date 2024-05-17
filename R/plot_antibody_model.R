
#' Antibody dependent boosting relationship
#'
#' Calculates the inferred antibody dependent boosting relationship from the MCMC chain
#' @param chain the MCMC chain
#' @param n number of samples to take
#' @param titres the vector of titres to calculate boosting values at
#' @return a data frame of quantiles for the inferred boost from different titre levels
#' @export
plot_antibody_dependent_boosting <- function(chain, n, titres = seq(0, 8, by = 0.1)) {
  samp_nos <- sample(unique(chain$samp_no), n)
  store <- matrix(nrow = n, ncol = length(titres))
  i <- 1
  for (samp in samp_nos) {
    pars <- as.numeric(chain[chain$samp_no == samp, ])
    names(pars) <- colnames(chain)
    mu <- pars["boost_long"] + pars["boost_short"]
    
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

#' Plot the antibody model
#' 
#' Plots the trajectory of the serosolver antibody model using specified parameters and optionally a specified antigenic map and infection history.
#' @inheritParams simulate_antibody_model
#' @param label_parameters if TRUE, labels the model parameters on the plot
#' @return a list with two ggplot objects, one showing the simulated antibody kinetics over time, stratified by biomarker ID, the other showing simulated antibody kinetics for each biomarker ID, stratified by sample time
#' @examples
#' plot_antibody_model(c("boost_long"=2,"boost_short"=3,"boost_delay"=1,"wane_short"=0.2,"wane_long"=0.01, "antigenic_seniority"=0,"cr_long"=0.1,"cr_short"=0.03), times=seq(1,25,by=1),infection_history=NULL,antigenic_map=example_antigenic_map)
#' @export
plot_antibody_model <- function(pars, 
                                times=NULL, 
                                infection_history=NULL, 
                                antigenic_map=NULL,
                                label_parameters=FALSE){
  y <- simulate_antibody_model(pars,times, infection_history,antigenic_map)
  y$biomarker_id_label <- paste0("Biomarker ID: ", y$biomarker_ids)
  y$sample_label <- paste0("Sample time: ", y$sample_times)
  
  y$biomarker_id_label <- factor(y$biomarker_id_label, levels=unique(paste0("Biomarker ID: ", y$biomarker_ids)))
  y$sample_label <- factor(y$sample_label, levels=unique(paste0("Sample time: ", y$sample_times)))
  
  if(is.null(infection_history)){
    infection_history <- 1
  }
  
  ## Add some labels and lines to show how the model parameters work
  parameter_labels_long <- data.frame(names=c("boost_long","boost_short","boost_delay"),
                                      biomarker_id_label = "Biomarker ID: 0",
                                      x_label = c(pars["boost_delay"] + 2, pars["boost_delay"] + 2,pars["boost_delay"]/2 + 1 ),
                                      y_label = c(pars["boost_long"]/2, (pars["boost_long"] + pars["boost_short"] + pars["boost_long"])/2,-0.1),
                                      
                                      xmin=c(pars["boost_delay"] + 1,pars["boost_delay"] + 1, 1),
                                      xmax=c(pars["boost_delay"] + 1,pars["boost_delay"] + 1,1 + pars["boost_delay"]),
                                      ymin=c(0,pars["boost_long"],0),
                                      ymax=c(pars["boost_long"],pars["boost_long"] + pars["boost_short"], 0)
                                      )
  
  p_long <- ggplot(y) + 
    geom_line(aes(x=sample_times,y=antibody_level)) + 
    geom_vline(data=data.frame(x=infection_history),linetype="dashed",aes(col="Infection",xintercept=x)) +
    facet_wrap(~biomarker_id_label) +
    scale_color_manual(name="",values=c("Infection"="grey40")) +
    scale_y_continuous(expand=c(0.03,0.03)) +
    scale_x_continuous(expand=c(0.03,0.03)) +
    xlab("Sample time") +
    ylab("Antibody level") +
    theme_pubr()+
    theme(legend.position="bottom")
  
  if(label_parameters){
    wane_short_x <- 1 + max(pars["boost_delay"] + (1-0.25/pars["boost_short"])/pars["wane_short"],0)
    wane_short_y <- pars["boost_short"]*0.25 + pars["boost_long"]
    wane_long_x <- wane_short_x + 5
    wane_long_y <- pars["boost_long"]*0.75
    
    p_long <- p_long + 
    geom_segment(data=parameter_labels_long,aes(x=xmin,xend=xmax,y=ymin,yend=ymax,group=names),linetype="dotted")+
    geom_text(data=parameter_labels_long, aes(label=names,x=x_label, y = y_label)) +
    geom_text(data=data.frame(label="wane_short", x=wane_short_x, y=wane_short_y,biomarker_id_label="Biomarker ID: 0"),aes(label=label,x=x,y=y)) +
    geom_text(data=data.frame(label="wane_long", x=wane_long_x, y=wane_long_y,biomarker_id_label="Biomarker ID: 0"),aes(label=label,x=x,y=y))
  }
  
  p_cr <- ggplot(y) + 
    geom_ribbon(aes(x=biomarker_ids,ymax=antibody_level,ymin=0),fill='grey70',col="black") + 
    geom_vline(data=data.frame(x=infection_history),linetype="dashed",aes(col="Infection",xintercept=x)) +
    facet_wrap(~sample_label) +
    scale_color_manual(name="", values=c("Infection"="grey40")) +
    scale_y_continuous(expand=c(0.03,0.03)) +
    scale_x_continuous(expand=c(0.03,0.03)) +
    xlab("Biomarker ID") +
    ylab("Antibody level") +
    theme_pubr()+
    theme(legend.position="bottom")
  
  return(list(p_long,p_cr))
}

#' Plots infection histories and antibody model 
#'
#' Given outputs from an MCMC run and the data used for fitting, generates an NxM matrix of plots where N is the number of individuals to be plotted and M is the range of sampling times. Where data are available, plots the observed antibody measurements and model predicted trajectories. Unlike plot_infection_histories_cross_sectional, places biomarker\_id on the x-axis and facets by sample time and individual.
#' @inheritParams get_antibody_level_predictions
#' @param known_infection_history nxm matrix of known infection histories
#' @param p_ncol integer giving the number of columns of subplots to create if using orientation = "longitudinal"
#' @param orientation either "cross-sectional" or "longitudinal"
#' @param subset_biomarker_ids if not NULL, then a vector giving the entries of biomarker_id to include in the longitudinal plot
#' @param settings if not NULL, list of serosolver settings as returned from the main serosolver function
#' @return a ggplot2 object
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_theta_chain)
#' data(example_inf_chain)
#' data(example_antibody_data)
#' data(example_antigenic_map)
#' data(example_par_tab)
#'
#' model_fit_plot <- plot_model_fits(example_theta_chain, example_inf_chain, example_antibody_data, 
#'                                            1:10, example_antigenic_map, example_par_tab,orientation="longitudinal")
#' }
#' @export
plot_model_fits <- function(chain, infection_histories, 
                            antibody_data=NULL,  
                            demographics=NULL,
                            individuals, 
                            par_tab=NULL,
                            antigenic_map=NULL, 
                            possible_exposure_times=NULL,
                            nsamp = 1000,
                            known_infection_history=NULL,
                            measurement_bias = NULL,
                            p_ncol=max(1,floor(length(individuals)/2)),
                            data_type=1,
                            expand_to_all_times=FALSE,
                            orientation="cross-sectional",
                            subset_biomarker_ids=NULL,
                            subset_biomarker_groups = NULL,
                            start_level="none",
                            settings=NULL
) {
  ## If the list of serosolver settings was included, use these rather than passing each one by one
  if(!is.null(settings)){
    message("Using provided serosolver settings list")
    if(is.null(antigenic_map)) antigenic_map <- settings$antigenic_map
    if(is.null(possible_exposure_times)) possible_exposure_times <- settings$possible_exposure_times
    if(is.null(measurement_bias)) measurement_bias <- settings$measurement_bias
    if(is.null(antibody_data)) antibody_data <- settings$antibody_data
    if(is.null(demographics)) demographics <- settings$demographics
    if(is.null(par_tab)) par_tab <- settings$par_tab
    if(is.null(start_level) | start_level == "none") start_level <- settings$start_level
    if(missing(data_type)) data_type <- settings$data_type
  }
  individuals <- individuals[order(individuals)]
  
  ## Setup antigenic map and exposure times
  ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
  if (!is.null(antigenic_map)) {
    possible_exposure_times_tmp <- unique(antigenic_map$inf_times) 
    if(is.null(possible_exposure_times)) {
      possible_exposure_times <- possible_exposure_times_tmp
    }
  } 
  if(class(start_level) == "character"){
  start_levels <- create_start_level_data(antibody_data %>% 
                                            dplyr::filter(individual %in% individuals),start_level,FALSE) %>% 
                                            dplyr::arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
  } else if(class(start_level) %in% c("tibble","data.frame")){
    start_levels <- start_levels
  } else {
    start_levels <- NULL
  }
  ## Generate antibody predictions
  antibody_preds <- get_antibody_level_predictions(
    chain, infection_histories, antibody_data, 
    demographics,
    individuals,
    antigenic_map, possible_exposure_times, 
    par_tab, nsamp, FALSE, 
    measurement_bias,
    expand_antibody_data=TRUE,
    expand_to_all_times=expand_to_all_times,
    data_type=data_type,
    start_level=start_levels
  )

  ## Use these antibody predictions and summary statistics on infection histories
  to_use <- antibody_preds$predicted_observations
  model_preds <- antibody_preds$predictions
  to_use$individual <- individuals[to_use$individual]

  inf_hist_densities <- antibody_preds$histories
  inf_hist_densities$xmin <- inf_hist_densities$variable-0.5
  inf_hist_densities$xmax <- inf_hist_densities$variable+0.5
  ## Subset infection history densities to not plot infections before sample time
  inf_hist_densities <- inf_hist_densities %>% 
    left_join(model_preds[model_preds$individual %in% individuals,c("individual","sample_time")] %>% dplyr::distinct(),by="individual",relationship="many-to-many") %>% 
    dplyr::filter(variable <= sample_time)
  
  if(is.null(par_tab)){
    measurement_ranges <- antibody_data %>% group_by(biomarker_group) %>% dplyr::summarize(min_measurement=min(measurement,na.rm=TRUE),
                                                                                           max_measurement=max(measurement,na.rm=TRUE))
  } else{
    if(!"biomarker_group" %in% colnames(par_tab)) par_tab$biomarker_group <- 1
    measurement_ranges <- par_tab %>% dplyr::filter(names %in% c("min_measurement","max_measurement")) %>% dplyr::select(names,values,biomarker_group) %>%
      pivot_wider(names_from=names,values_from=values)
  }
  max_x <- max(inf_hist_densities$variable) + 5
  time_range <- range(inf_hist_densities$variable)
  ## If provided, add true infection histories
  if(!is.null(known_infection_history)){
    known_infection_history <- known_infection_history[individuals,]
    rownames(known_infection_history) <- 1:nrow(known_infection_history)
    known_infection_history <- reshape2::melt(known_infection_history)
    colnames(known_infection_history) <- c("individual","variable","inf")
    known_infection_history <- known_infection_history[known_infection_history$inf == 1,]
    known_infection_history$variable <- possible_exposure_times[as.numeric(as.factor(known_infection_history$variable))]
    known_infection_history$individual <- individuals[known_infection_history$individual]
    expand_samples <- expand_grid(individual=individuals,sample_time=unique(to_use$sample_time))
    known_infection_history <- known_infection_history %>% left_join(expand_samples,by="individual",relationship="many-to-many") %>% filter(variable <= sample_time)
  }
  
  if(is.null(subset_biomarker_groups)){
    subset_biomarker_groups_use <- 1
  } else {
    subset_biomarker_groups_use <- subset_biomarker_groups
  }
  titre_pred_p <- NULL
  for(biomarker_group_use in subset_biomarker_groups_use){
    if(orientation=="cross-sectional"){
      if(!is.null(subset_biomarker_ids)){
        time_range <- range(subset_biomarker_ids)
      }
      
      p_tmp <- ggplot(to_use %>% dplyr::filter(biomarker_group == biomarker_group_use)) +
        geom_rect(data=inf_hist_densities%>% dplyr::cross_join(measurement_ranges)%>% dplyr::filter(biomarker_group == biomarker_group_use),
                  aes(xmin=xmin,xmax=xmax,fill=value,ymin=min_measurement-1,ymax=max_measurement+1))+
        geom_ribbon(aes(x=biomarker_id,ymin=lower, ymax=upper),alpha=0.4, fill="#009E73",linewidth=0.2)+
        geom_ribbon(data=model_preds[model_preds$individual %in% individuals,]%>% dplyr::filter(biomarker_group == biomarker_group_use), 
                    aes(x=biomarker_id,ymin=lower,ymax=upper),alpha=0.7,fill="#009E73",linewidth=0.2) + 
        geom_line(data=model_preds%>% dplyr::filter(biomarker_group == biomarker_group_use), aes(x=biomarker_id, y=median),linewidth=0.75,color="#009E73")+
        geom_rect(data=measurement_ranges,aes(ymin=max_measurement,ymax=max_measurement+1),xmin=0,xmax=max_x,fill="grey70")+
        geom_rect(data=measurement_ranges, aes(ymin=min_measurement-1,ymax=min_measurement),xmin=0,xmax=max_x,fill="grey70")
      
      
      if(!is.null(known_infection_history)){
        p_tmp <- p_tmp + geom_vline(data=known_infection_history,aes(xintercept=variable,linetype="Known infection")) +
          scale_linetype_manual(name="",values=c("Known infection"="dashed"))
      }
      min_measurement <- measurement_ranges %>% dplyr::filter(biomarker_group == biomarker_group_use) %>% dplyr::pull(min_measurement)
      max_measurement <- measurement_ranges %>% dplyr::filter(biomarker_group == biomarker_group_use) %>% dplyr::pull(max_measurement)
      breaks <- seq(floor(min_measurement), floor(max_measurement),by=2)
      
      p_tmp <- p_tmp +
        scale_x_continuous(expand=c(0.01,0.01)) +
        scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
        guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                                    barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
        geom_point(data=antibody_data %>% dplyr::filter(individual %in% individuals) %>%
                     dplyr::filter(biomarker_group == biomarker_group_use), aes(x=biomarker_id, y=measurement),shape=23, 
                   col="black",size=1)+
        ylab("log antibody level") +
        xlab("Time of antigen circulation") +
        theme_pubr()+
        theme(legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.margin = margin(-1,-1,-3,-1),
              axis.title=element_text(size=10),
              axis.text.x=element_text(angle=45,hjust=1,size=8),
              axis.text.y=element_text(size=8),
              plot.margin=margin(r=15,t=5,l=5))+
        coord_cartesian(xlim=time_range,ylim=c(min(breaks)-1, max(breaks)+1)) +
        scale_y_continuous(expand=c(0,0),breaks=breaks) +
        facet_grid(individual~sample_time)
    } else {
      
      if(!is.null(subset_biomarker_ids)){
        to_use <- to_use %>% dplyr::filter(biomarker_id %in% subset_biomarker_ids)
        model_preds <- model_preds %>% dplyr::filter(biomarker_id %in% subset_biomarker_ids)
        antibody_data <- antibody_data %>% dplyr::filter(biomarker_id %in% subset_biomarker_ids)
      }
        
      to_use$biomarker_id <- as.factor(to_use$biomarker_id)
      model_preds$biomarker_id <- as.factor(to_use$biomarker_id)
      antibody_data$biomarker_id <- as.factor(antibody_data$biomarker_id)
      p_tmp <- ggplot(to_use[to_use$individual %in% individuals,]%>% dplyr::filter(biomarker_group == biomarker_group_use)) +
        geom_rect(data=inf_hist_densities %>% dplyr::cross_join(measurement_ranges)%>% dplyr::filter(biomarker_group == biomarker_group_use) %>% select(-sample_time) %>% distinct(),
                  aes(xmin=xmin,xmax=xmax,alpha=value,ymin=min_measurement-1,ymax=max_measurement+1),fill="orange")+
        geom_ribbon(aes(x=sample_time,ymin=lower, ymax=upper,fill=biomarker_id,group=biomarker_id),alpha=0.1, linewidth=0.2)+
        geom_ribbon(data=model_preds[model_preds$individual %in% individuals,]%>% dplyr::filter(biomarker_group == biomarker_group_use), 
                    aes(x=sample_time,ymin=lower,ymax=upper,fill=biomarker_id,group=biomarker_id),alpha=0.25,linewidth=0.2) + 
        geom_line(data=model_preds%>% dplyr::filter(biomarker_group == biomarker_group_use), aes(x=sample_time, y=median,color=biomarker_id,group=biomarker_id),linewidth=0.75)+
        geom_rect(data=measurement_ranges%>% dplyr::filter(biomarker_group == biomarker_group_use),aes(ymin=max_measurement,ymax=max_measurement+1),xmin=0,xmax=max_x,fill="grey70")+
        geom_rect(data=measurement_ranges%>% dplyr::filter(biomarker_group == biomarker_group_use),aes(ymin=min_measurement-1,ymax=min_measurement),xmin=0,xmax=max_x,fill="grey70")
      
      
      if(!is.null(known_infection_history)){
        p_tmp <- p_tmp + geom_vline(data=known_infection_history,aes(xintercept=variable,linetype="Known infection")) +
          scale_linetype_manual(name="",values=c("Known infection"="dashed"))
      }
      
      min_measurement <- measurement_ranges %>% dplyr::filter(biomarker_group == biomarker_group_use) %>% dplyr::pull(min_measurement)
      max_measurement <- measurement_ranges %>% dplyr::filter(biomarker_group == biomarker_group_use) %>% dplyr::pull(max_measurement)
      breaks <- seq(floor(min_measurement), floor(max_measurement),by=2)
      
      p_tmp <- p_tmp +
        scale_x_continuous(expand=c(0.01,0.01)) +
        scale_alpha_continuous(range=c(0,1),name="Posterior probability of infection")+
        scale_fill_viridis_d(name="Biomarker ID") +
        scale_color_viridis_d(name="Biomarker ID") +
        geom_point(data=antibody_data %>% dplyr::filter(individual %in% individuals) %>%
                     dplyr::filter(biomarker_group == biomarker_group_use), aes(x=sample_time, y=measurement,col=biomarker_id),shape=23, 
                   size=1)+
        ylab("log antibody level") +
        xlab("Time of antigen circulation") +
        theme_pubr()+
        theme(legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.margin = margin(-1,-1,-3,-1),
              axis.title=element_text(size=10),
              axis.text.x=element_text(angle=45,hjust=1,size=8),
              axis.text.y=element_text(size=8),
              plot.margin=margin(r=15,t=5,l=5),
              legend.position="bottom")+
        coord_cartesian(xlim=range(to_use$sample_time),ylim=c(min(breaks)-1, max(breaks)+1)) +
        scale_y_continuous(expand=c(0,0),breaks=breaks) +
        facet_wrap(~individual,ncol=p_ncol)
    }
    titre_pred_p[[biomarker_group_use]] <- p_tmp
  }
  if(is.null(subset_biomarker_groups)){
    titre_pred_p <- titre_pred_p[[1]]
  }
  titre_pred_p
}

## From ggpubr
theme_pubr <- function (base_size = 12, base_family = "", border = FALSE, margin = TRUE, 
                          legend = c("top", "bottom", "left", "right", "none"), x.text.angle = 0) 
{
  half_line <- base_size/2
  if (!is.numeric(legend)) 
    legend <- match.arg(legend)
  if (x.text.angle > 5) 
    xhjust <- 1
  else xhjust <- NULL
  if (border) {
    panel.border <- element_rect(fill = NA, colour = "black", 
                                 size = 0.7)
    axis.line <- element_blank()
  }
  else {
    panel.border <- element_blank()
    axis.line = element_line(colour = "black", size = 0.5)
  }
  if (margin) 
    plot.margin <- margin(half_line, half_line, half_line, 
                          half_line)
  else plot.margin <- unit(c(0.5, 0.3, 0.3, 0.3), "mm")
  .theme <- theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = panel.border, panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = axis.line, 
          axis.text = element_text(color = "black"), legend.key = element_blank(), 
          strip.background = element_rect(fill = "#F2F2F2", 
                                          colour = "black", size = 0.7), plot.margin = plot.margin, 
          legend.position = legend, complete = TRUE)
  if (x.text.angle != 0) 
    .theme <- .theme + theme(axis.text.x = element_text(angle = x.text.angle, 
                                                        hjust = xhjust))
  .theme
}


#' Plots model predicted titers against observations
#'
#' @inheritParams plot_model_fits
#' @return a list with: 
#' \itemize{
#' \item a data frame with all posterior estimates for each observation; 
#' \item the proportion of observations captured by the 95% prediction intervals; 
#' \item a histogram comparing posterior median estimates to the observed data (note, this can be misleading for continuous data due to the zero-inflated observation model); 
#' \item a histogram comparing random posterior draws to the observed data (can be more reliable than posterior medians);
#' \item comparison of observations and all posterior medians and 95% prediction intervals
#' }
#' @family infection_history_plots
#' @export
plot_antibody_predictions <- function(chain, infection_histories, 
                                      antibody_data=NULL, 
                                      demographics=NULL,
                                      par_tab=NULL,
                                      antigenic_map=NULL, 
                                      possible_exposure_times=NULL,
                                      nsamp = 1000,
                                      measurement_bias = NULL,
                                      data_type=1,
                                      start_level="none",
                                      settings=NULL){
  ## If the list of serosolver settings was included, use these rather than passing each one by one
  if(!is.null(settings)){
    message("Using provided serosolver settings list")
    if(is.null(antigenic_map)) antigenic_map <- settings$antigenic_map
    if(is.null(measurement_bias)) measurement_bias <- settings$measurement_bias
    if(is.null(antibody_data)) antibody_data <- settings$antibody_data
    if(is.null(demographics)) demographics <- settings$demographics
    if(is.null(par_tab)) par_tab <- settings$par_tab
    if(is.null(start_level) | start_level == "none") start_level <- settings$start_level
    if(missing(data_type)) data_type <- settings$data_type
  }
  
  ## Setup antigenic map and exposure times
  ## Check if an antigenic map is provided. If not, then create a dummy map where all pathogens have the same position on the map
  if (!is.null(antigenic_map)) {
    possible_exposure_times_tmp <- unique(antigenic_map$inf_times) 
    if(is.null(possible_exposure_times)) {
      possible_exposure_times <- possible_exposure_times_tmp
    }
  } 
  if(class(start_level) == "character"){
    start_levels <- create_start_level_data(antibody_data,start_level,FALSE) %>% 
    dplyr::arrange(individual, biomarker_group, sample_time, biomarker_id, repeat_number)
  } else if(class(start_level) %in% c("tibble","data.frame")){
    start_levels <- start_levels
  } else {
    start_levels <- NULL
  }
  ## Generate antibody predictions
  x <- get_antibody_level_predictions(
    chain, infection_histories, antibody_data, individuals=unique(antibody_data$individual),
    antigenic_map, possible_exposure_times, 
    par_tab, nsamp, FALSE, 
    measurement_bias,
    expand_antibody_data=FALSE,
    expand_to_all_times=FALSE,
    data_type=data_type,
    start_level=start_levels,
    demographics=demographics,
    for_regression = TRUE
  )
  
  ## Use these antibody predictions and summary statistics on infection histories
  ## Find proportion of measurements within the 95% prediction intervals
  prop_correct <- x$predicted_observations %>% 
    mutate(correct = measurement >= lower & measurement <= upper) %>% group_by(correct) %>% tally() %>%
    pivot_wider(names_from=correct,values_from=n)
  if(!("FALSE" %in% colnames(prop_correct))) {
    prop_correct$`FALSE` <- 0
  }
  if(!("TRUE" %in% colnames(prop_correct))) {
    prop_correct$`TRUE` <- 0
  }
  prop_correct <- prop_correct %>%
    mutate(prop_correct1 = `TRUE`/(`FALSE` + `TRUE`)) %>% 
    pull(prop_correct1)
  
  
  ## Plot data against 95% prediction intervals and posterior medians
  p_compare_pointrange <- x$predicted_observations %>% arrange(measurement) %>%mutate(i=1:n()) %>% ggplot() + 
    geom_pointrange(aes(x=i,y=median,ymin=lower,ymax=upper,col="Predicted observation"),linewidth=0.25,size=0.25,alpha=0.5) + 
    geom_point(aes(x=i,y=measurement,col="Observation"),size=0.25) + 
    theme_minimal() + 
    scale_color_manual(name="",values=c("Predicted observation"="blue","Observation"="black")) +
    ylab("Antibody level") + 
    xlab("Measurement ID\n(ascending measurement value)") +
    theme(legend.position="bottom") +
    ggtitle("Comparison of observations to posterior median and 95% prediction intervals")
  
  ## Histograms of observations against predictions
  ## Posterior median
  p_hist_median <- ggplot(x$predicted_observations %>% mutate(i=1:n()) %>%
                            dplyr::select(i,measurement,median) %>%
                            rename(`Posterior median`=median,`Observation`=measurement) %>%
                            pivot_longer(-i)) +
    geom_histogram(aes(x=value,fill=name),position="dodge")+
    theme_minimal() + 
    scale_fill_manual(name="",values=c("Posterior median"="blue","Observation"="black")) +
    theme(legend.position="bottom") +
    xlab("Antibody level") +
    ylab("Count") +
    ggtitle("Distribution of measurements vs. posterior median predictions")
  
  ## Some random draws
  rand_draws <- sample(1:ncol(x$all_predictions_obs),min(9, ncol(x$all_predictions_obs)))
  tmp_pred_obs <- x$all_predictions_obs[,rand_draws]
  colnames(tmp_pred_obs) <- rand_draws
  p_hist_draws <- x$predicted_observations %>% select(measurement) %>% bind_cols(tmp_pred_obs) %>%
    mutate(i=1:n()) %>%
    pivot_longer(-c(i,measurement)) %>%
    mutate(name=paste0("Posterior draw: ", name)) %>%
    ggplot() + 
    geom_histogram(aes(x=measurement,fill="Observation"),alpha=0.5,col="grey10",
                   linewidth=0.1) + 
    geom_histogram(aes(x=value,fill="Predicted"),alpha=0.5,col="grey10",linewidth=0.1) +
    scale_fill_manual(name="",values=c("Predicted"="blue","Observation"="black")) +
    
    facet_wrap(~name) +
    theme_minimal() +
    theme(legend.position="bottom") +
    xlab("Antibody level") +
    ylab("Count") +
    ggtitle("Predicted observations vs. measurements from random posterior draws")
  
  
  return(list("all_predictions"=x$predicted_observations,
              "proportion_correct"=paste0("Proportion of observations within 95% prediction intervals: ", signif(prop_correct,5)),
              "p_hist_median"=p_hist_median,
              "p_hist_draws"=p_hist_draws,
              "p_pointrange"=p_compare_pointrange))
  
}

#' Plots estimated antibody kinetics model
#'
#' @inheritParams plot_model_fits
#' @param solve_times vector of times to solve model over
#' @return a ggplot2 object giving model-predicted antibody level and predicted observations over time since infection
#' @family infection_history_plots
#' @export
plot_estimated_antibody_model <- function(chain, 
                                          antibody_data=NULL, 
                                          demographics = NULL,
                                          antigenic_map=NULL,
                                          possible_exposure_times=NULL, 
                                          par_tab=NULL,
                                          nsamp = 1000, 
                                          measurement_bias = NULL,
                                          solve_times = seq(1,30,by=1),
                                          data_type=1,
                                          settings=NULL,
                                          by_group=TRUE){
  ## If the list of serosolver settings was included, use these rather than passing each one by one
  if(!is.null(settings)){
    message("Using provided serosolver settings list")
    if(is.null(antigenic_map)) antigenic_map <- settings$antigenic_map
    if(is.null(possible_exposure_times)) possible_exposure_times <- settings$possible_exposure_times
    if(is.null(measurement_bias)) measurement_bias <- settings$measurement_bias
    if(is.null(antibody_data)) antibody_data <- settings$antibody_data
    if(is.null(demographics)) demographics <- settings$demographics
    if(is.null(par_tab)) par_tab <- settings$par_tab
    if(missing(data_type)) data_type <- settings$data_type
  }
  
  par_tab <- add_scale_pars(par_tab,antibody_data,demographics)
  
  ## Get unique demographic groups from full data set, not just the subset
  if(!is.null(demographics)){
    demographic_groups <- create_demographic_table(demographics,par_tab)
  } else {
    demographic_groups <- create_demographic_table(antibody_data,par_tab)
  }
  
  ## Need to align the iterations of the two MCMC chains
  ## and choose some random samples
  
  ## Convert samp_no and chainno to a single samp_no index
  if(!("chain_no" %in% colnames(chain))){
    chain$chain_no <- 1
  }
  
  chain <- chain %>% dplyr::group_by(chain_no,samp_no) %>% 
    dplyr::mutate(samp_no = cur_group_id()) %>% dplyr::ungroup() %>% 
    dplyr::mutate(chain_no = 1) %>% arrange(samp_no)
  
  samps <- unique(chain$samp_no)
  nsamp <- min(nsamp, length(unique(chain$samp_no)))
  
  labels <- rep("", nrow(demographic_groups))
  for(i in 1:nrow(demographic_groups)) 
    for(j in 1:ncol(demographic_groups)) 
      labels[i] <- paste0(labels[i], paste0(colnames(demographic_groups)[j], ":", demographic_groups[i,j],";"))
  
  demographic_groups_plot <- demographic_groups %>% mutate(individual = 1:n())
  demographic_groups_plot$Group <- labels
  
  ## Create fake antibody data for all individuals and times
  full_antibody_data <- antibody_data %>% 
    select(biomarker_group,biomarker_id) %>% distinct() %>% 
    cross_join(demographic_groups_plot) %>% 
    expand_grid(sample_time=solve_times) %>% 
    mutate(birth = 1, repeat_number=1,measurement=0)
  n_indiv <- length(unique(full_antibody_data$individual))
  unique_biomarker_groups <- unique(full_antibody_data$biomarker_group)
  
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
  
  tmp_samp <- sample(samps, nsamp)
  ## See the function in posteriors.R
  model_func <- create_posterior_func(par_tab, full_antibody_data, antigenic_map, possible_exposure_times,
                                      prior_version=2,
                                      measurement_bias = measurement_bias, function_type = 4,
                                      antibody_level_before_infection=FALSE,
                                      data_type=data_type,start_level="none",
                                      demographics=demographics,demographic_groups=demographic_groups
  )
  
  predicted_titres <- observed_predicted_titres <- matrix(nrow = nrow(full_antibody_data), ncol = nsamp)
  samp_record <- numeric(nsamp)
  
  inf_hist_test <- matrix(0, nrow=n_indiv, ncol = length(possible_exposure_times))
  inf_hist_test[,1] <- 1
  ## For each sample, take values for theta and infection histories and simulate titres
  for (i in 1:nsamp) {
    index <- tmp_samp[i]
    pars <- get_index_pars(chain, samp_no=index,chain_no=1)
    pars <- pars[!(names(pars) %in% c("posterior_prob", "likelihood", "prior_prob",
                                      "samp_no", "total_infections", "chain_no"
    ))]
    names(pars) <- par_tab$names
    
    predicted_titres[, i] <- model_func(pars, inf_hist_test)
    for(biomarker_group in unique_biomarker_groups){
      observed_predicted_titres[which(full_antibody_data$biomarker_group == biomarker_group),i] <- add_noise(predicted_titres[which(full_antibody_data$biomarker_group == biomarker_group),i], pars, NULL, NULL,data_type=data_type[biomarker_group])
    }
    samp_record[i] <- index
  }
  
  colnames(predicted_titres) <- tmp_samp
  
  ## Get 95% credible interval and means
  dat2 <- t(apply(predicted_titres, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))))
  
  ## Get 95% credible interval and means of observations
  obs_dat <- t(apply(observed_predicted_titres, 1, function(x) c(mean(x), quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))))
  
  dat2 <- as.data.frame(dat2)
  obs_dat <- as.data.frame(obs_dat)
  
  colnames(dat2) <- colnames(obs_dat) <- c("mean","lower", "lower_50", "median", "upper_50", "upper")
  dat2 <- cbind(full_antibody_data, dat2)
  obs_dat <- cbind(full_antibody_data, obs_dat)
  
  ## Get blanks for ranges
  ranges <- par_tab[par_tab$names %in% c("min_measurement","max_measurement"),c("names","values","biomarker_group")] %>% pivot_wider(names_from=names,values_from=values)  

  
  if(by_group){
    dat2$biomarker_id <- as.factor(dat2$biomarker_id)
    obs_dat$biomarker_id <- as.factor(obs_dat$biomarker_id)
    
    p1 <- ggplot(dat2) + 
      geom_hline(data=ranges, aes(yintercept=min_measurement),linetype="dotted",linewidth=0.25) +
      geom_hline(data=ranges, aes(yintercept=max_measurement),linetype="dotted",linewidth=0.25) +
      geom_ribbon(data=obs_dat, aes(x=sample_time,ymin=lower,ymax=upper,y=mean,fill=biomarker_id,group=biomarker_id),alpha=0.1) +
      geom_ribbon(aes(x=sample_time,ymin=lower,ymax=upper,y=mean,fill=biomarker_id,group=biomarker_id),alpha=0.5) +
      geom_line(aes(x=sample_time,y=mean,group=biomarker_id,col=biomarker_id)) + 
      scale_y_continuous(expand=c(0.03,0.03)) +
      scale_x_continuous(expand=c(0.03,0.03),limits=range(solve_times)) +
      scale_fill_viridis_d(name="Biomarker ID: ") +
      scale_color_viridis_d(name="Biomarker ID: ") +
      xlab("Time since infection") +
      ylab("Antibody level") +
      theme_pubr()+
      theme(legend.position="bottom") +
      facet_wrap(paste0("Biomarker group: ", biomarker_group)~Group)
  } else {
    p1 <- ggplot(dat2) + 
      geom_hline(data=ranges, aes(yintercept=min_measurement),linetype="dotted",linewidth=0.25) +
      geom_hline(data=ranges, aes(yintercept=max_measurement),linetype="dotted",linewidth=0.25) +
      geom_ribbon(data=obs_dat, aes(x=sample_time,ymin=lower,ymax=upper,y=mean,fill=Group,group=Group),alpha=0.1) +
      geom_ribbon(aes(x=sample_time,ymin=lower,ymax=upper,y=mean,fill=Group,group=Group),alpha=0.5) +
      geom_line(aes(x=sample_time,y=mean,group=Group,col=Group)) + 
      scale_y_continuous(expand=c(0.03,0.03)) +
      scale_x_continuous(expand=c(0.03,0.03),limits=range(solve_times)) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() +
      xlab("Time since infection") +
      ylab("Antibody level") +
      theme_pubr()+
      theme(legend.position="bottom") +
      facet_grid(paste0("Biomarker ID: ", biomarker_id)~paste0("Biomarker group: ", biomarker_group))
  }
  
  
  return(p1)
}

