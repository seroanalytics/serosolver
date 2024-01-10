#' Antibody dependent boosting relationship
#'
#' Calculates the inferred antibody dependent boosting relationship from the MCMC chain
#' @param chain the MCMC chain
#' @param n number of samples to take
#' @param titres the vector of titres to calculate boosting values at
#' @return a data frame of quantiles for the inferred boost from different titre levels
#' @export
plot_antibody_dependent_boosting <- function(chain, n, titres = seq(0, 8, by = 0.1)) {
  sampnos <- sample(unique(chain$sampno), n)
  store <- matrix(nrow = n, ncol = length(titres))
  i <- 1
  for (samp in sampnos) {
    pars <- as.numeric(chain[chain$sampno == samp, ])
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
#' @param p_ncol integer giving the number of columns of subplots to create if using orientation = "cross-sectional"
#' @param orientation either "longitudinal" or "cross-sectional"
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
                                          individuals, antigenic_map=NULL, 
                                          possible_exposure_times=NULL, par_tab,
                                          nsamp = 100,
                                          measurement_indices_by_time = NULL,
                            p_ncol=length(individuals)/2,data_type=1,
                            orientation="longitudinal") {
  individuals <- individuals[order(individuals)]
  ## Generate antibody predictions
  titre_preds <- get_antibody_level_predictions(
    chain, infection_histories, antibody_data, individuals,
    antigenic_map, possible_exposure_times, 
    par_tab, nsamp, FALSE, 
    measurement_indices_by_time,
    expand_antibody_data=TRUE,
    data_type=data_type
  )
  
  ## Use these antibody predictions and summary statistics on infection histories
  to_use <- titre_preds$predicted_observations
  model_preds <- titre_preds$predictions
  to_use$individual <- individuals[to_use$individual]
  
  inf_hist_densities <- titre_preds$histories
  inf_hist_densities$xmin <- inf_hist_densities$variable-0.5
  inf_hist_densities$xmax <- inf_hist_densities$variable+0.5
  
  max_measurement <- max(antibody_data$measurement,na.rm=TRUE)
  min_measurement <- min(antibody_data$measurement,na.rm=TRUE)
  
  max_x <- max(inf_hist_densities$variable) + 5
  time_range <- range(inf_hist_densities$variable)
  
  if(orientation=="longitudinal"){
    titre_pred_p <- ggplot(to_use) +
      geom_rect(data=inf_hist_densities,
                aes(xmin=xmin,xmax=xmax,fill=value),ymin=min_measurement-1,ymax=max_measurement+2)+
      geom_ribbon(aes(x=biomarker_id,ymin=lower, ymax=upper),alpha=0.4, fill="#009E73",size=0.2)+
      geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
                  aes(x=biomarker_id,ymin=lower,ymax=upper),alpha=0.7,fill="#009E73",size=0.2) + 
      geom_line(data=model_preds, aes(x=biomarker_id, y=median),linewidth=0.75,color="#009E73")+
      geom_rect(ymin=max_measurement,ymax=max_measurement+2,xmin=0,xmax=max_x,fill="grey70")+
      geom_rect(ymin=min_measurement-2,ymax=min_measurement,xmin=0,xmax=max_x,fill="grey70")+
      scale_x_continuous(expand=c(0,0)) +
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
    titre_pred_p <- ggplot(to_use) +
      geom_rect(data=inf_hist_densities,
                aes(xmin=xmin,xmax=xmax,fill=value),ymin=min_measurement-1,ymax=max_measurement+2)+
      geom_ribbon(aes(x=sample_time,ymin=lower, ymax=upper),alpha=0.25, fill="#009E73",size=0.2)+
      geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
                  aes(x=sample_time,ymin=lower,ymax=upper),alpha=0.5,fill="#009E73",size=0.2) + 
      geom_line(data=model_preds, aes(x=sample_time, y=median),linewidth=0.75,color="#009E73")+
      geom_rect(ymin=max_measurement,ymax=max_measurement+2,xmin=0,xmax=max_x,fill="grey70")+
      geom_rect(ymin=min_measurement-2,ymax=min_measurement,xmin=0,xmax=max_x,fill="grey70")+
      scale_x_continuous(expand=c(0,0)) +
      scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
      guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                                  barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
      geom_point(data=antibody_data[antibody_data$individual %in% individuals,], aes(x=sample_time, y=measurement),shape=23, 
                 col="black",size=1,fill=viridis(1)[1])+
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
      coord_cartesian(ylim=c(min_measurement,max_measurement+1),xlim=range(possible_exposure_times)) +
      scale_y_continuous(breaks=seq(min_measurement,max_measurement+2,by=2)) +
      facet_wrap(~individual,ncol=p_ncol)
  }
  titre_pred_p
}
