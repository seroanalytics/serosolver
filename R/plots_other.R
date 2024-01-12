#' Useful plot for looking at simulated data
#'
#' Plots measured titres and known infection histories for all individuals, faceted by sample time (multi-antigen panel) or biomarker_id variable (longitudinal single antigen)
#' @param antibody_data the data frame of antibody data
#' @param infection_histories the infection history matrix
#' @param possible_exposure_times the vector of times at which individuals could be infected
#' @param n_indivs how many individuals to plot
#' @param start_inf if not NULL, plots the infection history matrix used as the starting point in the MCMC chain
#' @param study_design default "multi-antigen" facets by sample time. "single-antigen" gives sample time on the x-axis and colours by biomarker_id
#' @return a ggplot object
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_antibody_data)
#' data(example_inf_hist)
#' data(example_antigenic_map)
#' possible_exposure_times <- example_antigenic_map$inf_times
#' plot_antibody_data(example_antibody_data, example_inf_hist, possible_exposure_times, 5)
#' }
#' @export
plot_antibody_data <- function(antibody_data, infection_histories, 
                      possible_exposure_times, 
                      n_indivs, start_inf = NULL,
                      study_design="multi-antigen"){
    indivs <- unique(antibody_data$individual)
    infection_history <- as.data.frame(cbind(indivs, infection_histories))
    colnames(infection_history) <- c("individual", possible_exposure_times)
    melted_inf_hist <- reshape2::melt(infection_history, id.vars = "individual")
    melted_inf_hist$variable <- as.numeric(as.character(melted_inf_hist$variable))
    melted_inf_hist <- melted_inf_hist[melted_inf_hist$value > 0, ]
    tmp <- unique(antibody_data[, c("individual", "sample_time")])
    melted_inf_hist <- merge(melted_inf_hist, tmp)
    melted_inf_hist <- melted_inf_hist[melted_inf_hist$variable <= melted_inf_hist$sample_time, ]
    samps <- sample(unique(antibody_data$individual), n_indivs)

    if("biomarker_group" %in% colnames(antibody_data)){
        antibody_data$biomarker_group <- as.factor(antibody_data$biomarker_group)
    } else {
        antibody_data$biomarker_group <- 1
    }
    
    if (study_design == "multi-antigen") {
        p1 <- ggplot(antibody_data[antibody_data$individual %in% samps, ]) +
            geom_point(aes(x = as.integer(biomarker_id), y = measurement, col=biomarker_group)) +
            geom_vline(data = melted_inf_hist[melted_inf_hist$individual %in% samps, ], aes(xintercept = variable,linetype="Known infection time"), col = "grey10") +
            theme_bw() +
            xlab("Antigen") +
            facet_grid(individual ~ sample_time)
    } else {
        p1 <- ggplot(antibody_data[antibody_data$individual %in% samps, ]) +
            geom_point(aes(x = sample_time, y = measurement, col=biomarker_id, col=biomarker_group)) +
            geom_vline(data = melted_inf_hist[melted_inf_hist$individual %in% samps, ], 
                       aes(xintercept = variable,linetype="Known infection time"), col = "grey10") +
            theme_bw() +
            xlab("Antigen circulation time") +
            facet_wrap(~individual)
    }
    
    if (!is.null(start_inf)) {
        start_inf_hist <- as.data.frame(cbind(indivs, start_inf))
        colnames(start_inf_hist) <- c("individual", possible_exposure_times)
        melted_start_hist <- reshape2::melt(start_inf_hist, id.vars = "individual")
        melted_start_hist$variable <- as.numeric(as.character(melted_start_hist$variable))
        melted_start_hist <- melted_start_hist[melted_start_hist$value > 0, ]
        p1 <- p1 + geom_vline(data = melted_start_hist[melted_start_hist$individual %in% samps, ], aes(xintercept = variable), col = "blue", linetype = "dashed")
    }
    p1 <- p1 + ylab("log antibody level") + scale_color_viridis_d() + scale_linetype_manual(name="",values=c("Known infection time"="dashed"))
    
    p1 <- p1 + 
      theme(legend.position="bottom",
            axis.text.x=element_text(size=6),
            axis.text.y=element_text(size=6),
            axis.title=element_text(size=8))
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
