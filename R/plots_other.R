#' Plot raw data
#'
#' Plots measured antibody measurements and known infection histories for all individuals, faceted by sample time (multi-antigen panel) or biomarker_id variable (longitudinal single antigen)
#' @param antibody_data the data frame of antibody data
#' @param possible_exposure_times the vector of times at which individuals could be infected
#' @param n_indivs integer of how many individuals to plot, or vector of which individuals to plot
#' @param infection_histories the infection history matrix
#' @param study_design default "multi-antigen" facets by sample time. "single-antigen" gives sample time on the x-axis and colours by biomarker_id
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
                      study_design="multi-antigen"){
    indivs <- unique(antibody_data$individual)
    
    if(length(n_indivs) == 1){
      samps <- sample(unique(antibody_data$individual), n_indivs)
    } else {
      samps <- n_indivs
    }

    max_measurement <- max(antibody_data$measurement)
    min_measurement <- min(antibody_data$measurement)
    
    max_x <- max(possible_exposure_times)
    time_range <- range(possible_exposure_times)
    
    if("biomarker_group" %in% colnames(antibody_data)){
        antibody_data$biomarker_group <- as.factor(antibody_data$biomarker_group)
    } else {
        antibody_data$biomarker_group <- 1
    }
    p1 <- ggplot(antibody_data[antibody_data$individual %in% samps, ]) +
      geom_rect(ymin=max_measurement,ymax=max_measurement+2,xmin=0,xmax=max_x,fill="grey70") +
      geom_rect(ymin=min_measurement-2,ymax=min_measurement,xmin=0,xmax=max_x,fill="grey70") 
      
    if (study_design == "multi-antigen") {
        p1 <- p1 + 
            geom_point(aes(x = as.integer(biomarker_id), y = measurement, col=biomarker_group),shape=23, 
                       col="black",size=1)+
            facet_grid(individual ~ sample_time)
    } else {
        p1 <- p1 + 
            geom_point(aes(x = sample_time, y = measurement, col=biomarker_id, col=biomarker_group),shape=23, 
                       col="black",size=1) +
            theme_bw() +
            facet_wrap(~individual)
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
      scale_x_continuous(expand=c(0,0)) +
      ylab("log antibody level") + 
      xlab("Time of antigen circulation") +
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
      coord_cartesian(ylim=c(min_measurement,max_measurement+1),xlim=time_range) +
      scale_y_continuous(breaks=seq(min_measurement,max_measurement+2,by=2)) +
      scale_color_viridis_d() + 
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
