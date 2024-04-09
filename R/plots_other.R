#' Plot raw data
#'
#' Plots measured antibody measurements and known infection histories for all individuals, faceted by sample time (multi-antigen panel) or biomarker_id variable (longitudinal single antigen)
#' @param antibody_data the data frame of antibody data
#' @param possible_exposure_times the vector of times at which individuals could be infected
#' @param n_indivs integer of how many individuals to plot, or vector of which individuals to plot
#' @param infection_histories the infection history matrix
#' @param study_design default "cross-sectional" facets by sample time. "longitudinal" gives sample time on the x-axis and colours by `biomarker_id`
#' @param measurement_ranges data frame or tibble stating for each `biomarker_group` the `min_measurement` and `max_measurement`. If NULL, this is extracted from `antibody_data`
#' @return a ggplot object
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_antibody_data)
#' data(example_inf_hist)
#' data(example_antigenic_map)
#' possible_exposure_times <- example_antigenic_map$inf_times
#' plot_antibody_data(example_antibody_data, possible_exposure_times, 5, example_inf_hist)
#' }
#' @export
plot_antibody_data <- function(antibody_data, 
                               possible_exposure_times, 
                               n_indivs,     
                               infection_histories=NULL,                   
                      study_design="cross-sectional",
                      measurement_ranges=NULL){
    indivs <- unique(antibody_data$individual)
    
    if(length(n_indivs) == 1){
      samps <- sample(unique(antibody_data$individual), n_indivs)
    } else {
      samps <- n_indivs
    }

    ## Note that this might not be the actual range
    if(is.null(measurement_ranges)){
      measurement_ranges <- antibody_data %>% group_by(biomarker_group) %>% dplyr::summarize(min_measurement=min(measurement,na.rm=TRUE),
                                                                                            max_measurement=max(measurement,na.rm=TRUE))
    }
    
    max_x <- max(possible_exposure_times)
    time_range <- range(possible_exposure_times)
    
    if("biomarker_group" %in% colnames(antibody_data)){
        antibody_data$biomarker_group <- as.factor(antibody_data$biomarker_group)
    } else {
        antibody_data$biomarker_group <- 1
    }
    p1 <- ggplot(antibody_data[antibody_data$individual %in% samps, ]) +
      geom_rect(data=measurement_ranges,aes(ymin=max_measurement,ymax=max_measurement+2),xmin=0,xmax=max_x,fill="grey70") +
      geom_rect(data=measurement_ranges,aes(ymin=min_measurement-2,ymax=min_measurement),xmin=0,xmax=max_x,fill="grey70") 
      
    if (study_design == "cross-sectional") {
        p1 <- p1 + 
            geom_point(aes(x = as.integer(biomarker_id), y = measurement, col=biomarker_group),shape=23, 
                       col="black",size=1)+
            facet_grid(individual ~ sample_time) +
        xlab("Time of antigen circulation")
    } else {
        p1 <- p1 + 
            geom_point(aes(x = sample_time, y = measurement, col=as.factor(biomarker_id)),size=1) +
            xlab("Sample time") +
            theme_bw() +
            facet_grid(biomarker_group~individual)
    }
    
    if(!is.null(infection_histories)) {
      infection_history <- as.data.frame(cbind(indivs, infection_histories))
      colnames(infection_history) <- c("individual", possible_exposure_times)
      melted_inf_hist <- reshape2::melt(infection_history, id.vars = "individual")
      melted_inf_hist$variable <- as.numeric(as.character(melted_inf_hist$variable))
      melted_inf_hist <- melted_inf_hist[melted_inf_hist$value > 0, ]
      tmp <- unique(antibody_data[, c("individual", "sample_time")])
      melted_inf_hist <- merge(melted_inf_hist, tmp)
      melted_inf_hist <- melted_inf_hist[melted_inf_hist$variable <= melted_inf_hist$sample_time, ]
      
      p1 <- p1 + geom_vline(data = melted_inf_hist[melted_inf_hist$individual %in% samps, ], 
                            aes(xintercept = variable,linetype="Known infection time"), 
                            col = "orange",linetype="dashed")
    }
    
    p1 <- p1 +
      scale_x_continuous(expand=c(0.05,0.05)) +
      ylab("log antibody level") + 
      theme_pubr()+
      theme(
        legend.position="bottom",
        legend.title=element_text(size=7),
            legend.text=element_text(size=7),
            legend.margin = margin(-1,-1,-3,-1),
            axis.title=element_text(size=10),
            axis.text.x=element_text(angle=45,hjust=1,size=8),
            axis.text.y=element_text(size=8),
            plot.margin=margin(r=15,t=5,l=5))+
      coord_cartesian(xlim=time_range) +
      scale_y_continuous(expand=c(0,0)) + 
      scale_color_viridis_d(name="Biomarker ID") + 
      scale_linetype_manual(name="",values=c("Known infection time"="dashed"))
    return(p1)
}

#' @family theta_plots
#' @export
plot_2d_density <- function(chain, par1, par2) {
    ggplot(chain) +
        stat_density_2d(aes_string(
            x = par1, y = par2,
            alpha = "stat(level)",
            col = "as.factor(chain_no)",
            fill = "as.factor(chain_no)"
        ),
        geom = "polygon", size = 0.2
        ) +
        scale_alpha_continuous(range = c(0.01, 0.3))
}

#' Plot time between serum samples
#'
#' @param antibody_data the data frame of antibody data, including labels for individuals and time sample was taken
#' @return a ggplot2 object
#' @family theta_plots
#' @examples
#' \dontrun{
#' data(example_antibody_data)
#' plot_samples_distances(example_antibody_data)
#' }
#' @export
plot_samples_distances <- function(antibody_data) {
    sample_times <- unique(antibody_data[, c("individual", "sample_time")])
    distances <- ddply(sample_times, ~individual, function(x) {
        if (nrow(x) < 2) {
            y <- 0
        } else {
            y <- diff(x$sample_time)
        }
        y
    })
    ggplot(distances) + geom_histogram(aes(x = V1), binwidth = 1) + theme_bw() + xlab("Time points between samples")
}

#' @export
plot_posteriors_theta <- function(chain,par_tab){
  par_tab_tmp <- par_tab[par_tab$fixed == 0,]
  par_indices <- c("samp_no","chain_no",par_tab_tmp[par_tab_tmp$par_type == 1,"names"],"total_infections","likelihood","prior_prob","posterior_prob")

  chain <- chain %>% as_tibble()
  chain$chain_no <- as.factor(chain$chain_no)
  
  p_trace_par <- chain[,par_indices] %>%
    pivot_longer(-c(samp_no, chain_no)) %>%
    ggplot() + 
    geom_line(aes(x=samp_no,y=value,col=chain_no)) + 
    facet_wrap(~name,scales="free_y") +
    theme_bw()
  
  p_density_par <- chain[,par_indices] %>%
    pivot_longer(-c(samp_no, chain_no)) %>%
    ggplot() + 
    geom_density(aes(x=value,fill=chain_no),alpha=0.25) +
    facet_wrap(~name,scales="free") +
    theme_bw()
  
  p_trace_phi <- p_density_phi <- NULL
  if("phi" %in% par_tab$names & "phi" %in% colnames(chain)){
    phi_indices <- 1:nrow(par_tab_tmp[par_tab_tmp$names == "phi",])
    phi_indices <- c("samp_no","chain_no",c("phi",paste0("phi.",phi_indices[1:(length(phi_indices)-1)])))
    p_trace_phi <- chain[,phi_indices] %>%
      pivot_longer(-c(samp_no, chain_no)) %>%
      ggplot() + 
      geom_line(aes(x=samp_no,y=value,col=chain_no)) + 
      facet_wrap(~name,scales="free_y") +
      theme_bw()
    p_density_phi <- chain[,phi_indices] %>%
      pivot_longer(-c(samp_no, chain_no)) %>%
      ggplot() + 
      geom_density(aes(x=value,fill=chain_no),alpha=0.25) +
      facet_wrap(~name,scales="free") +
      theme_bw()
  }
  p_trace_rho <- p_density_rho <- NULL
  if("rho" %in% par_tab$names& "rho" %in% colnames(chain)){
    rho_indices <- 1:nrow(par_tab_tmp[par_tab_tmp$names == "rho",])
    rho_indices <- c("samp_no","chain_no",c("rho",paste0("rho.",rho_indices[1:(length(rho_indices)-1)])))
    p_trace_rho <- chain[,rho_indices] %>%
      pivot_longer(-c(samp_no, chain_no)) %>%
      ggplot() + 
      geom_line(aes(x=samp_no,y=value,col=chain_no)) + 
      facet_wrap(~name,scales="free_y") +
      theme_bw()
    p_density_rho <- chain[,rho_indices] %>%
      pivot_longer(-c(samp_no, chain_no)) %>%
      ggplot() + 
      geom_density(aes(x=value,fill=chain_no),alpha=0.25) +
      facet_wrap(~name,scales="free") +
      theme_bw()
  }
    
  list(p_trace_par, p_density_par,
       p_trace_phi,p_density_phi,
       p_trace_rho,p_density_rho)
}

#' @export
plot_mcmc_diagnostics <- function(location, par_tab, burnin){
  chains <- load_mcmc_chains(location=location,par_tab=par_tab,burnin=burnin,estimated_only=TRUE)
  
  par_medians <- sapply(chains$theta_chain[,!(colnames(chains$theta_chain) %in% c("samp_no","chain_no","posterior_prob","likelihood","prior_prob"))], function(x) median(x))
  par_means <- sapply(chains$theta_chain[,!(colnames(chains$theta_chain) %in% c("samp_no","chain_no","posterior_prob","likelihood","prior_prob"))], function(x) mean(x))
  par_lower95 <- sapply(chains$theta_chain[,!(colnames(chains$theta_chain) %in% c("samp_no","chain_no","posterior_prob","likelihood","prior_prob"))], function(x) quantile(x,0.025))
  par_upper95 <- sapply(chains$theta_chain[,!(colnames(chains$theta_chain) %in% c("samp_no","chain_no","posterior_prob","likelihood","prior_prob"))], function(x) quantile(x,0.975))
  
  miss_cols <- colnames(chains$theta_chain)
  miss_cols <- miss_cols[!(miss_cols %in% c("samp_no","chain_no","likelihood","prior_prob","posterior_prob"))]
  chains1 <- lapply(chains$theta_list_chains, function(x){
    x <- x %>% select(all_of(miss_cols)) %>% as.mcmc()
  })

  ess <- effectiveSize(chains1)
  
  par_estimates <- data.table(names=names(par_medians),median=par_medians,
                              mean=par_means,lower95_CrI=par_lower95,upper95_CrI=par_upper95,
                              ess=ess)
  
  if(length(chains$theta_list_chains) > 1){
    gelman_res <- gelman.diag(chains1)
    par_estimates <- cbind(par_estimates, gelman_res$psrf)
    colnames(par_estimates)[7:8] <- c("Rhat point estimate","Rhat upper CI")
  } else {
    gelman_res <- "Cannot calculate Rhat with only 1 chain"
  }
  
  
  
  p_inf_hists <- plot_posteriors_infhist(chains$inf_chain,example_antigenic_map$inf_times,n_alive=NULL)
  p_thetas <- plot_posteriors_theta(chains$theta_chain,par_tab)
  
  list(
    theta_estimates=par_estimates,
    p_thetas=p_thetas,
    inf_hist_estimates=p_inf_hists$estimates,
    p_inf_hists=p_inf_hists[which(names(p_inf_hists) != "estimates")]
      )
}
