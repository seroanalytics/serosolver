
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


#' Plots infection histories and antibody model 
#'
#' Given outputs from an MCMC run and the data used for fitting, generates an NxM matrix of plots where N is the number of individuals to be plotted and M is the range of sampling times. Where data are available, plots the observed antibody measurements and model predicted trajectories. Unlike plot_infection_histories_cross_sectional, places biomarker_id on the x-axis and facets by sample time and individual.
#' @inheritParams get_antibody_level_predictions
#' @param known_infection_history nxm matrix of known infection histories
#' @param p_ncol integer giving the number of columns of subplots to create if using orientation = "longitudinal"
#' @param orientation either "cross-sectional" or "longitudinal"
#' @param subset_biomarker_ids if not NULL, then a vector giving the entries of biomarker_id to include in the longitudinal plot
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
plot_model_fits <- function(chain, infection_histories, antibody_data,
                                          individuals, par_tab,
                            antigenic_map=NULL, 
                                          possible_exposure_times=NULL,
                                          nsamp = 1000,
  known_infection_history=NULL,
                                          measurement_indices_by_time = NULL,
                            p_ncol=length(individuals)/2,
                            data_type=1,
                            expand_to_all_times=FALSE,
                            orientation="cross-sectional",
                            subset_biomarker_ids=NULL) {
  individuals <- individuals[order(individuals)]
  if(is.null(possible_exposure_times)){
    possible_exposure_times <- unique(antigenic_map$inf_times)
 }
  ## Generate antibody predictions
  antibody_preds <- get_antibody_level_predictions(
    chain, infection_histories, antibody_data, individuals,
    antigenic_map, possible_exposure_times, 
    par_tab, nsamp, FALSE, 
    measurement_indices_by_time,
    expand_antibody_data=TRUE,
    expand_to_all_times=expand_to_all_times,
    data_type=data_type
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
    left_join(model_preds[model_preds$individual %in% individuals,c("individual","sample_time")] %>% dplyr::distinct(),by="individual") %>% 
    dplyr::filter(variable <= sample_time)
  
  max_measurement <- max(antibody_data$measurement,na.rm=TRUE)
  min_measurement <- min(antibody_data$measurement,na.rm=TRUE)
  
  max_x <- max(inf_hist_densities$variable) + 5
  time_range <- range(inf_hist_densities$variable)
  
  ## If provided, add true infection histories
  if(!is.null(known_infection_history)){
    known_infection_history <- known_infection_history[individuals,]
    known_infection_history <- reshape2::melt(known_infection_history)
    colnames(known_infection_history) <- c("individual","variable","inf")
    known_infection_history <- known_infection_history[known_infection_history$inf == 1,]
    known_infection_history$variable <- possible_exposure_times[known_infection_history$variable]
    known_infection_history$individual <- individuals[known_infection_history$individual]
    expand_samples <- expand_grid(individual=individuals,sample_time=unique(to_use$sample_time))
    known_infection_history <- known_infection_history %>% left_join(expand_samples,by="individual") %>% filter(variable <= sample_time)
  }
  
  if(orientation=="cross-sectional"){
    titre_pred_p <- ggplot(to_use) +
      geom_rect(data=inf_hist_densities,
                aes(xmin=xmin,xmax=xmax,fill=value),ymin=min_measurement-1,ymax=max_measurement+2)+
      geom_ribbon(aes(x=biomarker_id,ymin=lower, ymax=upper),alpha=0.4, fill="#009E73",size=0.2)+
      geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
                  aes(x=biomarker_id,ymin=lower,ymax=upper),alpha=0.7,fill="#009E73",size=0.2) + 
      geom_line(data=model_preds, aes(x=biomarker_id, y=median),linewidth=0.75,color="#009E73")+
      geom_rect(ymin=max_measurement,ymax=max_measurement+2,xmin=0,xmax=max_x,fill="grey70")+
      geom_rect(ymin=min_measurement-2,ymax=min_measurement,xmin=0,xmax=max_x,fill="grey70")
    
    
    if(!is.null(known_infection_history)){
      titre_pred_p <- titre_pred_p + geom_vline(data=known_infection_history,aes(xintercept=variable,linetype="Known infection")) +
        scale_linetype_manual(name="",values=c("Known infection"="dashed"))
    }
    
    titre_pred_p <- titre_pred_p +
      scale_x_continuous(expand=c(0.01,0.01)) +
      scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
      guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                                  barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
      geom_point(data=antibody_data[antibody_data$individual %in% individuals,], aes(x=biomarker_id, y=measurement),shape=23, 
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
      coord_cartesian(ylim=c(min_measurement,max_measurement+1),xlim=time_range) +
      scale_y_continuous(breaks=seq(min_measurement,max_measurement+2,by=2)) +
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
    
      
    titre_pred_p <- ggplot(to_use[to_use$individual %in% individuals,]) +
      geom_ribbon(aes(x=sample_time,ymin=lower, ymax=upper,fill=biomarker_id,group=biomarker_id),alpha=0.25, size=0.2)+
      #geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
      #            aes(x=sample_time,ymin=lower,ymax=upper,fill=biomarker_id,group=biomarker_id),alpha=0.5,size=0.2) + 
      geom_line(data=model_preds, aes(x=sample_time, y=median,color=biomarker_id,group=biomarker_id),linewidth=0.75)+
      geom_rect(ymin=max_measurement,ymax=max_measurement+2,xmin=0,xmax=max_x,fill="grey70")+
      geom_rect(ymin=min_measurement-2,ymax=min_measurement,xmin=0,xmax=max_x,fill="grey70")
    
    
    if(!is.null(known_infection_history)){
      titre_pred_p <- titre_pred_p + geom_vline(data=known_infection_history,aes(xintercept=variable,linetype="Known infection")) +
        scale_linetype_manual(name="",values=c("Known infection"="dashed"))
    }
    
    titre_pred_p <- titre_pred_p +
      scale_x_continuous(expand=c(0.01,0.01)) +
      #scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
     # guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
      #                            barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
      scale_fill_viridis_d(name="Biomarker ID") +
      scale_color_viridis_d(name="Biomarker ID") +
      geom_point(data=antibody_data[antibody_data$individual %in% individuals,], aes(x=sample_time, y=measurement,col=biomarker_id),shape=23, 
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
      coord_cartesian(ylim=c(min_measurement,max_measurement+1),xlim=range(to_use$sample_time)) +
      scale_y_continuous(breaks=seq(min_measurement,max_measurement+2,by=2)) +
      facet_wrap(~individual,ncol=p_ncol)
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
