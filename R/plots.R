#' Formatted quantiles
#'
#' Given a vector of MCMC samples, generates and formats the desired quantile estimates
#' @param x the vector to summarise
#' @param sig_f how many significant figures to print
#' @param qs the vector of quantiles
#' @param as_text if TRUE, formats nicely as text rather than a vector of numbers
#' @return the formatted quantiles
#' @examples
#' data(example_theta_chain)
#' x <- example_theta_chain$mu
#' generate_quantiles(x)
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
#' @param antigenic_map (optional) a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. See \code{\link{example_antigenic_map}}
#' @param strain_isolation_times (optional) if no antigenic map is specified, this argument gives the vector of times at which individuals can be infected
#' @param par_tab the table controlling the parameters in the MCMC chain
#' @param nsamp number of samples to take from posterior
#' @param add_residuals if true, returns an extra output summarising residuals between the model prediction and data
#' @param mu_indices vector of integers. for random effects on boosting parameter, mu. If random mus are included in the parameter table, this vector specifies which mu to use for each circulation year. For example, if years 1970-1976 have unique boosting, then mu_indices should be c(1,2,3,4,5,6). If every 3 year block shares has a unique boosting parameter, then this should be c(1,1,1,2,2,2)
#' @param measurement_indices_by_time default NULL, optional vector giving the index of `measurement_bias` that each strain uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clusters
#' @param for_res_plot TRUE/FALSE value. If using the output of this for plotting of residuals, returns the actual data points rather than summary statistics
#' @param expand_titredat TRUE/FALSE value. If TRUE, solves titre predictions for all possible infection times. If left FALSE, then only solves for the infections times at which a titre against the circulating virus was measured in titre_dat.
#' @param titre_before_infection TRUE/FALSE value. If TRUE, solves titre predictions, but gives the predicted titre at a given time point BEFORE any infection during that time occurs.
#' @param data_type integer, currently accepting 1 or 2. Set to 1 for discretized, bounded data, or 2 for continuous, bounded data. Note that with 2, MIN_TITRE must be set.
#' @return a list with the titre predictions (95% credible intervals, median and multivariate posterior mode) and the probabilities of infection for each individual in each epoch
#' @examples
#' \dontrun{
#' data(example_theta_chain)
#' data(example_inf_chain)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' data(example_par_tab)
#'
#' y <- get_titre_predictions(example_theta_chain, example_inf_chain, example_titre_dat,
#'                           unique(example_titre_dat$individual), example_antigenic_map,
#'                           example_par_tab,expand_titredat = FALSE)
#' }
#' @export
get_titre_predictions <- function(chain, infection_histories, titre_dat,
                                  individuals, antigenic_map=NULL,
                                  strain_isolation_times=NULL, par_tab,
                                  nsamp = 100, add_residuals = FALSE,
                                  mu_indices = NULL,
                                  measurement_indices_by_time = NULL,
                                  for_res_plot = FALSE, expand_titredat = FALSE,
                                  titre_before_infection=FALSE, titres_for_regression=FALSE,
                                  data_type=1){
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
    if (!is.null(antigenic_map)) {
        strain_isolation_times <- unique(antigenic_map$inf_times) # How many strains are we testing against and what time did they circulate
    } else {
        antigenic_map <- data.frame("x_coord"=1,"y_coord"=1,"inf_times"=strain_isolation_times)
    }
    nstrain <- length(strain_isolation_times)
    n_indiv <- length(individuals)

    ## Empty data structures to save output to
    infection_history_dens <- NULL
    tmp_samp <- sample(samps, nsamp)
    
    ## See the function in posteriors.R
    titre_dat1 <- titre_dat
    
    if (expand_titredat) {
        titre_dat1 <- expand.grid(
            individual = unique(titre_dat$individual),
            samples = unique(titre_dat$samples),
            obs_type=unique(titre_dat$obs_type),
            titre = 0, run = 1
        )
        titre_dat2 <- unique(titre_dat[, c("individual", "virus", "group", "DOB")])
        titre_dat1 <- merge(titre_dat1, titre_dat2)
        titre_dat1 <- titre_dat1[
            order(titre_dat1$group, titre_dat1$individual, titre_dat1$samples, titre_dat1$virus),
            c("individual", "samples","obs_type", "virus", "titre", "run", "group", "DOB")
        ]
    }
    model_func <- create_posterior_func(par_tab, titre_dat1, antigenic_map, 100,
                                        mu_indices = mu_indices,version=2,
                                        measurement_indices_by_time = measurement_indices_by_time, function_type = 4,
                                        titre_before_infection=titre_before_infection,
                                        data_type=data_type
                                        )

    predicted_titres <- residuals <- residuals_floor <- 
        observed_predicted_titres <- matrix(nrow = nrow(titre_dat1), ncol = nsamp)
    samp_record <- numeric(nsamp)


    ## For each sample, take values for theta and infection histories and simulate titres
    inf_hist_all <- list(nsamp)
    for (i in 1:nsamp) {
        index <- tmp_samp[i]
        pars <- get_index_pars(chain, index)
        pars <- pars[!(names(pars) %in% c("lnlike", "likelihood", "prior_prob",
                                            "sampno", "total_infections", "chain_no"
                                        ))]
                                        ## pars <- pars[names(pars) %in% par_tab$names]
        tmp_inf_hist <- infection_histories[infection_histories$sampno == index, ]
        tmp_inf_hist <- as.matrix(Matrix::sparseMatrix(i = tmp_inf_hist$i, j = tmp_inf_hist$j, x = tmp_inf_hist$x, dims = c(n_indiv, nstrain)))
        predicted_titres[, i] <- model_func(pars, tmp_inf_hist)
        observed_predicted_titres[,i] <- add_noise(predicted_titres[,i], pars, NULL, NULL,data_type=data_type)
        inf_hist_all[[i]] <- tmp_inf_hist
        ## Get residuals between observations and predictions
        residuals[, i] <- titre_dat1$titre - floor(predicted_titres[, i])
        residuals_floor[,i] <- titre_dat1$titre - observed_predicted_titres[,i]
        samp_record[i] <- index
    }
    colnames(predicted_titres) <- tmp_samp

    ## If generating for residual plot, can return now
    if (for_res_plot) {
        return(list(residuals, samp_record, titre_dat1,
                    predicted_titres,
                    observed_predicted_titres,
                    residuals_floor))
    }

                                        #residuals <- cbind(titre_dat1, residuals)

    ## Get 95% credible interval and means
    dat2 <- t(apply(predicted_titres, 1, function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))))
    
    ## Get 95% credible interval and means of observations
    obs_dat <- t(apply(observed_predicted_titres, 1, function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))))
    
    residuals <- t(apply(residuals, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))
    residuals <- cbind(titre_dat1, residuals)

    ## Find multivariate posterior mode estimate from the chain
    best_pars <- get_best_pars(chain)
    best_pars <- best_pars[!(names(best_pars) %in% c(
                                                       "lnlike", "likelihood", "prior_prob",
                                                       "sampno", "total_infections", "chain_no"
                                                   ))]
                                        #best_pars <- best_pars[names(best_pars) %in% par_tab$names]
    best_I <- chain$sampno[which.max(chain$lnlike)]
    best_inf <- infection_histories[infection_histories$sampno == best_I, ]
    best_inf <- as.matrix(Matrix::sparseMatrix(i = best_inf$i, j = best_inf$j, x = best_inf$x, dims = c(n_indiv, nstrain)))

    ## Generate trajectory for best parameters
    best_traj <- model_func(best_pars, best_inf)
    best_residuals <- titre_dat1$titre - floor(best_traj)
    best_residuals <- cbind(titre_dat1, best_residuals, "sampno" = best_I)
    dat2 <- as.data.frame(dat2)
    obs_dat <- as.data.frame(obs_dat)
    
    colnames(dat2) <- colnames(obs_dat) <- c("lower", "lower_50", "median", "upper_50", "upper")
    dat2$max <- best_traj
    dat2 <- cbind(titre_dat1, dat2)
    obs_dat <- cbind(titre_dat1, obs_dat)
    tmp_inf_chain <- data.table(subset(infection_histories, sampno %in% tmp_samp))

    ## Get infection history density for each individual and each epoch
    data.table::setkey(tmp_inf_chain, "i", "j")
    infection_history_dens <- tmp_inf_chain[, list(V1 = sum(x) / length(tmp_samp)), by = key(tmp_inf_chain)]
    infection_history_dens$j <- strain_isolation_times[infection_history_dens$j]
    colnames(infection_history_dens) <- c("individual", "variable", "value")
    infection_history_final <- infection_history_dens
    best_inf <- data.frame(best_inf)
    best_inf$individual <- 1:nrow(best_inf)
    best_inf$individual <- individuals[best_inf$individual]

    dat2$individual <- individuals[dat2$individual]
    infection_history_final$individual <- individuals[infection_history_final$individual]
    if(titres_for_regression){
        return(list("all_predictions"=predicted_titres, "all_inf_hist"=inf_hist_all,
                    "summary_titres"=dat2,"best_inf_hist"=best_inf, "predicted_observations"=obs_dat)) 
    }
    
    if (add_residuals) {
        result <- list("predictions" = dat2, "histories" = infection_history_final, 
                       "residuals" = residuals, "bestRes" = best_residuals,"best_infhist"=best_inf,
                       "predicted_observations"=obs_dat)
    } else {
        result <- list("predictions" = dat2, "histories" = infection_history_final,
                       "best_infhist"=best_inf, "predicted_observations"=obs_dat)
    }
    return(result)
}

#' Plots infection histories and titre model fits longitudinal
#'
#' Given outputs from an MCMC run and the data used for fitting, generates an NxM matrix of plots where N is the number of individuals to be plotted and M is the range of sampling times. Where data are available, plots the observed titres and model predicted trajectories. Unlike plot_infection_histories, places virus on the x-axis and facets by sample and individual.
#' @inheritParams get_titre_predictions
#' @return a ggplot2 object
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_theta_chain)
#' data(example_inf_chain)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' data(example_par_tab)
#'
#' model_fit_plot <- plot_infection_histories_long(example_theta_chain, example_inf_chain, example_titre_dat, 
#'                                            1:10, example_antigenic_map, example_par_tab)
#' }
#' @export
plot_infection_histories_long <- function(chain, infection_histories, titre_dat,
                                     individuals, antigenic_map=NULL, 
                                     strain_isolation_times=NULL, par_tab,
                                     nsamp = 100,
                                     mu_indices = NULL,
                                     measurement_indices_by_time = NULL,
                                     data_type=1) {
    individuals <- individuals[order(individuals)]
    ## Generate titre predictions
    titre_preds <- get_titre_predictions(
        chain, infection_histories, titre_dat, individuals,
        antigenic_map, strain_isolation_times, 
        par_tab, nsamp, FALSE, mu_indices,
        measurement_indices_by_time,
        expand_titredat=TRUE,
        data_type=data_type
    )
    
    ## Use these titre predictions and summary statistics on infection histories
    to_use <- titre_preds$predicted_observations
    model_preds <- titre_preds$predictions
    to_use$individual <- individuals[to_use$individual]
    
    inf_hist_densities <- titre_preds$histories
    inf_hist_densities$xmin <- inf_hist_densities$variable-0.5
    inf_hist_densities$xmax <- inf_hist_densities$variable+0.5
    
    max_titre <- max(titre_dat$titre,na.rm=TRUE)
    min_titre <- min(titre_dat$titre,na.rm=TRUE)
    
    max_x <- max(inf_hist_densities$variable) + 5
    time_range <- range(inf_hist_densities$variable)
    titre_pred_p <- ggplot(to_use) +
        geom_rect(data=inf_hist_densities,
                  aes(xmin=xmin,xmax=xmax,fill=value),ymin=min_titre-1,ymax=max_titre+2)+
        geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),alpha=0.4, fill="#009E73",size=0.2)+
        geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
                    aes(x=virus,ymin=lower,ymax=upper),alpha=0.7,fill="#009E73",size=0.2) + 
        geom_line(data=model_preds, aes(x=virus, y=median),linewidth=0.75,color="#009E73")+
        geom_rect(ymin=max_titre,ymax=max_titre+2,xmin=0,xmax=max_x,fill="grey70")+
        geom_rect(ymin=min_titre-2,ymax=min_titre,xmin=0,xmax=max_x,fill="grey70")+
        scale_x_continuous(expand=c(0,0)) +
        scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
        guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                                    barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
        geom_point(data=titre_dat[titre_dat$individual %in% individuals,], aes(x=virus, y=titre),shape=23, 
                   col="black",size=1)+
        ylab("log titre") +
        xlab("Time of virus circulation") +
        theme_pubr()+
        theme(legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.margin = margin(-1,-1,-3,-1),
              axis.title=element_text(size=10),
              axis.text.x=element_text(angle=45,hjust=1,size=8),
              axis.text.y=element_text(size=8),
              plot.margin=margin(r=15,t=5,l=5))+
        coord_cartesian(ylim=c(min_titre,max_titre+1),xlim=time_range) +
        scale_y_continuous(breaks=seq(min_titre,max_titre+2,by=2)) +
        facet_grid(individual~samples)
    titre_pred_p
}

#' Plots infection histories and titre model fits
#'
#' Given outputs from an MCMC run and the data used for fitting, generates an NxM matrix of plots where N is the number of individuals to be plotted and M is the range of sampling times. Where data are available, plots the observed titres and model predicted trajectories
#' @inheritParams get_titre_predictions
#' @return a ggplot2 object
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_theta_chain)
#' data(example_inf_chain)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' data(example_par_tab)
#'
#' model_fit_plot <- plot_infection_histories(example_theta_chain, example_inf_chain, example_titre_dat, 
#'                                            1:10, example_antigenic_map, example_par_tab)
#' }
#' @export
plot_infection_histories <- function(chain, infection_histories, titre_dat,
                                     individuals, antigenic_map=NULL, 
                                     strain_isolation_times=NULL, par_tab,
                                     nsamp = 100,
                                     mu_indices = NULL,
                                     measurement_indices_by_time = NULL,
                                     p_ncol=length(individuals)/2,
                                     data_type=1) {
    individuals <- individuals[order(individuals)]
    ## Generate titre predictions
    titre_preds <- get_titre_predictions(
        chain, infection_histories, titre_dat, individuals,
        antigenic_map, strain_isolation_times, 
        par_tab, nsamp, FALSE, mu_indices,
        measurement_indices_by_time,
        expand_titredat=TRUE,data_type=data_type
    )

    ## Use these titre predictions and summary statistics on infection histories
    to_use <- titre_preds$predicted_observations
    model_preds <- titre_preds$predictions
        to_use$individual <- individuals[to_use$individual]
    
    inf_hist_densities <- titre_preds$histories
    inf_hist_densities$xmin <- inf_hist_densities$variable-0.5
    inf_hist_densities$xmax <- inf_hist_densities$variable+0.5
    
    max_titre <- max(titre_dat$titre,na.rm=TRUE)
    min_titre <- min(titre_dat$titre,na.rm=TRUE)
    
    max_x <- max(inf_hist_densities$variable) + 5
    time_range <- range(inf_hist_densities$variable)
    
    titre_pred_p <- ggplot(to_use) +
        geom_rect(data=inf_hist_densities,
                  aes(xmin=xmin,xmax=xmax,fill=value),ymin=min_titre-1,ymax=max_titre+2)+
        geom_ribbon(aes(x=samples,ymin=lower, ymax=upper),alpha=0.25, fill="#009E73",size=0.2)+
        geom_ribbon(data=model_preds[model_preds$individual %in% individuals,], 
                    aes(x=samples,ymin=lower,ymax=upper),alpha=0.5,fill="#009E73",size=0.2) + 
        geom_line(data=model_preds, aes(x=samples, y=median),linewidth=0.75,color="#009E73")+
        geom_rect(ymin=max_titre,ymax=max_titre+2,xmin=0,xmax=max_x,fill="grey70")+
        geom_rect(ymin=min_titre-2,ymax=min_titre,xmin=0,xmax=max_x,fill="grey70")+
        scale_x_continuous(expand=c(0,0)) +
        scale_fill_gradient(low="white",high="#D55E00",limits=c(0,1),name="Posterior probability of infection")+
        guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,label.position = "bottom",
                                    barwidth=10,barheight = 0.5, frame.colour="black",ticks=FALSE)) +
        geom_point(data=titre_dat[titre_dat$individual %in% individuals,], aes(x=samples, y=titre),shape=23, 
                   col="black",size=1,fill=viridis(1)[1])+
        ylab("log titre") +
        xlab("Time of virus circulation") +
        theme_pubr()+
        theme(legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.margin = margin(-1,-1,-3,-1),
              axis.title=element_text(size=10),
              axis.text.x=element_text(angle=45,hjust=1,size=8),
              axis.text.y=element_text(size=8),
              plot.margin=margin(r=15,t=5,l=5))+
        coord_cartesian(ylim=c(min_titre,max_titre+1),xlim=range(strain_isolation_times)) +
        scale_y_continuous(breaks=seq(min_titre,max_titre+2,by=2)) +
        facet_wrap(~individual,ncol=p_ncol)
    titre_pred_p
}

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
#' ## Load in exaple data
#' data(example_inf_chain)
#' data(example_antigenic_map)
#' data(example_titre_dat)
#'
#' strain_isolation_times <- example_antigenic_map$inf_times
#' ## Setup known attack rates
#' n_alive <- get_n_alive(example_titre_dat, strain_isolation_times)
#' n_infs <- colSums(example_inf_hist)
#' known_ar <- n_infs/n_alive
#' known_ar <- data.frame("j"=strain_isolation_times,"AR"=known_ar,"group"=1)
#' 
#' ## Setup known infection histories
#' known_inf_hist <- data.frame(example_inf_hist)
#' colnames(known_inf_hist) <- strain_isolation_times
#' 
#' n_alive_group <- get_n_alive_group(example_titre_dat, strain_isolation_times,melt_dat = TRUE)
#' n_alive_group$j <- strain_isolation_times[n_alive_group$j]
#' all_plots <- plot_posteriors_infhist(example_inf_chain, strain_isolation_times, n_alive_group,
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
    if (is.null(inf_chain$group)) {
        inf_chain$group <- 1
    }
    if(samples > length(unique(inf_chain$sampno))) {
        stop("Error - number of samples requested is greater than number of MCMC samples provided")
    }

    ## Pads the chain with 0 entries (ie. non-infection events)
    if (pad_chain) inf_chain <- pad_inf_chain(inf_chain)

    ## Thin the chain a bit to increase speed for plots
    thin_chain <- inf_chain[inf_chain$sampno %in% sample(unique(inf_chain$sampno), samples, replace = T), ]

    ## Trace plot by time
    time_plot <- plot_infection_history_chains_time(thin_chain, 0, sample(1:length(years), 10), n_alive=NULL, FALSE)
    ## Trace plot by indiv
    indiv_plot <- plot_infection_history_chains_indiv(thin_chain, 0, 1:10, FALSE)
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



#' Plot inferred posteriors theta
#'
#' Produces and saves estimated posterior distributions for the antibody kinetics parameters
#' @param chain the full MCMC chain to generate titre trajectories from
#' @param par_tab the table controlling the parameters in the MCMC chain
#' @param burnin if not already discarded, discard burn in from chain (takes rows where sampno > burnin)
#' @param samples how many samples from the chain to take
#' @param calculate_ess if TRUE, calculates the ESS for all free parameters
#' @param plot_corr if TRUE, returns a pairwise correlation plot of free parameters
#' @param save_plots if TRUE, directly saves the plots as svgs
#' @param plot_mcmc if TRUE, plots the MCMC chain traces
#' @param save_loc the full directory path of where to save plots
#' @return a list of ggplot objects and a data frame of estimates
#' @family theta_plots
#' @examples 
#' \dontrun{
#' data(example_theta_chain)
#' data(example_par_tab)
#' plot_posteriors_theta(example_theta_chain,example_par_tab,samples=100)
#' }
#' @export
plot_posteriors_theta <- function(chain,
                                  par_tab,
                                  burnin = 0,
                                  samples = 100,
                                  calculate_ess = TRUE,
                                  plot_corr = TRUE,
                                  save_plots = FALSE,
                                  plot_mcmc = TRUE,
                                  save_loc = "") {
    if (is.null(chain$chain_no)) {
        chain$chain_no <- 1
    }
    if(samples > length(unique(chain$sampno))) {
        stop("Error - number of samples requested is greater than number of MCMC samples provided")
    }
    ## Combined chain
    free_chain <- chain[chain$sampno > burnin, ]

    ## Get quantiles and create table of results
    ## These are from Adam's summaries...
    free_chain$sigma1drop <- free_chain$mu * free_chain$sigma1
    free_chain$sigma2drop <- free_chain$mu_short * free_chain$sigma2
    free_chain$wane_titre <- free_chain$mu_short * free_chain$wane
    free_chain$errorCorrect1 <- pnorm(3, mean = 1.5, sd = free_chain$error) - pnorm(2, mean = 2.5, sd = free_chain$error)
    free_chain$errorCorrect2 <- pnorm(4, mean = 1.5, sd = free_chain$error) - pnorm(1, mean = 2.5, sd = free_chain$error)


    thin_free_chain <- free_chain[sample(1:nrow(free_chain), samples, replace = TRUE), ]

    parameter <- c(
        par_tab[which(par_tab$fixed == 0), "names"],
        "sigma1drop", "sigma2drop", "wane_titre", "errorCorrect1", "errorCorrect2"
    )


    parameter <- sapply(unique(parameter), function(x) {
        n_pars <- length(which(parameter == x))
        if (n_pars > 1) {
            to_append <- c("", paste0(".", seq(1, n_pars - 1)))
            parameter[which(parameter == x)] <- paste0(parameter[which(parameter == x)], to_append)
        } else {
            parameter[which(parameter == x)]
        }
    })
    parameter <- unlist(parameter)
    names(parameter) <- NULL
    parameter <- intersect(parameter, colnames(free_chain))

    par_names <- par_tab$names
    parameter_table_names <- sapply(unique(par_names), function(x) {
        n_pars <- length(which(par_names == x))
        if (n_pars > 1) {
            to_append <- c("", paste0(".", seq(1, n_pars - 1)))
            par_names[which(par_names == x)] <- paste0(par_names[which(par_names == x)], to_append)
        } else {
            par_names[which(par_names == x)]
        }
    })
    par_tab$names <- unlist(parameter_table_names)
    par_tab_free <- par_tab[par_tab$fixed == 0, c("names", "values")]
    colnames(par_tab_free)[2] <- "par_tab_value"

    thin_free_chain <- thin_free_chain[, parameter]
    free_chain <- free_chain[, parameter]

    results <- data.frame("estimate" = apply(thin_free_chain, 2, function(x) generate_quantiles(x)))
    results$names <- parameter
    results_unformatted <- as.data.frame(t(apply(thin_free_chain, 2, function(x) generate_quantiles(x, as_text = FALSE))))
    results_unformatted$names <- parameter
    results_unformatted <- merge(par_tab_free, results_unformatted, by = c("names"), all = TRUE)
    results_unformatted$correct <- results_unformatted$par_tab_value < results_unformatted$`97.5%` & results_unformatted$par_tab_value > results_unformatted$`2.5%`

    all_results <- merge(results, results_unformatted, by = "names")


    ## Plot correlations
    corP <- NULL
    if (plot_corr) {
        corP <- GGally::ggpairs(thin_free_chain) + theme_bw()

        if (save_plots) {
            to.svg(print(corP), paste0(save_loc, "corPlot.svg"))
            to.svg(coda::autocorr.plot(thin_free_chain), paste0(save_loc, "autocorr.svg"))
                                        # to.pdf(print(corP), paste0(save_loc,"corPlot.pdf"))
                                        # to.pdf(coda::autocorr.plot(free_chain),paste0(save_loc,"autocorr.pdf"))
        }
    }


    ## Calculate ESS for all free parameters
    all_ess <- NULL
    if (calculate_ess) {
        all_ess <- coda::effectiveSize(free_chain)
        all_ess <- data.frame("names" = names(all_ess), ESS = all_ess)
        all_results <- merge(all_results, all_ess)
    }

    densP <- traceP <- NULL
    if (plot_mcmc) {
        densP <- bayesplot::mcmc_hist(thin_free_chain) + theme()
        traceP <- bayesplot::mcmc_trace(thin_free_chain) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
        if (save_plots) {
            to.svg(print(densP), paste0(save_loc, "densities.svg"))
            to.svg(print(traceP), paste0(save_loc, "traces.svg"))
                                        # to.pdf(print(densP), paste0(save_loc,"densities.pdf"))
                                        # to.pdf(print(traceP), paste0(save_loc,"traces.pdf"))
        }
    }
    if (save_plots) write.table(all_results, paste0(save_loc, "estimates.csv"), row.names = FALSE, sep = ",")
    return(list(results = all_results, corPlot = corP, densP = densP, traceP = traceP))
}

#' Get posterior information infection histories
#'
#' Finds the median, mean and 95% credible intervals for the attack rates and total number of infections per individual
#' @param solve_cumulative if TRUE, finds the cumulative infection histories for each individual. This takes a while, so is left FALSE by default.
#' @inheritParams plot_posteriors_infhist
#' @return a list of data frames with summary statistics
#' @family infection_history_plots
#' @examples
#' data(example_inf_chain)
#' data(example_antigenic_map)
#' data(example_titre_dat)
#' data(example_inf_hist)
#' strain_isolation_times <- example_antigenic_map$inf_times
#' ## Find number alive in each time period
#' n_alive <- get_n_alive(example_titre_dat, strain_isolation_times)
#' ## Get actual number of infections per time
#' n_infs <- colSums(example_inf_hist)
#' ## Create data frame of true ARs
#' known_ar <- n_infs/n_alive
#' known_ar <- data.frame("j"=strain_isolation_times,"AR"=known_ar,"group"=1)
#' ## Get true infection histories
#' known_inf_hist <- data.frame(example_inf_hist)
#' colnames(known_inf_hist) <- strain_isolation_times
#'
#' ## Need to get group specific n_alive and adjust to correct time frame 
#' n_alive_group <- get_n_alive_group(example_titre_dat, strain_isolation_times,melt_dat = TRUE)
#' n_alive_group$j <- strain_isolation_times[n_alive_group$j]
#' results <- calculate_infection_history_statistics(example_inf_chain, 0, strain_isolation_times,
#'                                                   n_alive=n_alive_group, known_ar=known_ar,
#'                                                   known_infection_history=known_inf_hist)
#' @export
calculate_infection_history_statistics <- function(inf_chain, burnin = 0, years = NULL,
                                                   n_alive = NULL, known_ar = NULL,
                                                   group_ids = NULL,
                                                   known_infection_history = NULL,
                                                   solve_cumulative=FALSE) {
    inf_chain <- inf_chain[inf_chain$sampno > burnin, ]
    if (is.null(inf_chain$chain_no)) {
        inf_chain$chain_no <- 1
    }
    message("Padding inf chain...\n")
    inf_chain <- pad_inf_chain(inf_chain)
    message("Done\n")

    if (!is.null(group_ids)) {
        inf_chain <- merge(inf_chain, data.table(group_ids))
    } else {
        inf_chain$group <- 1
    }

    message("Calculating by time summaries...\n")
    data.table::setkey(inf_chain, "group", "j", "sampno", "chain_no")
    n_inf_chain <- inf_chain[, list(total_infs = sum(x)), by = key(inf_chain)]


    if (!is.null(years)) {
        n_inf_chain$j <- years[n_inf_chain$j]
    }
    
    if (!is.null(n_alive)) {
        n_inf_chain <- merge(n_inf_chain, n_alive, by = c("j", "group"))
        n_inf_chain$total_infs <- n_inf_chain$total_infs / n_inf_chain$n_alive
        n_inf_chain[is.nan(n_inf_chain$total_infs), "total_infs"] <- 0
    }
    data.table::setkey(n_inf_chain, "group", "sampno", "chain_no")
    n_inf_chain[, cumu_infs := cumsum(total_infs), by = key(n_inf_chain)]
    gelman_res_j <- ddply(n_inf_chain, .(group,j), function(tmp_chain){
        tmp_chain_mcmc <- split(as.data.table(tmp_chain), by=c("chain_no"))
        tmp_chain_mcmc <- lapply(tmp_chain_mcmc, function(x) as.mcmc(x[,c("total_infs")]))
        tmp_chain_mcmc <- as.mcmc.list(tmp_chain_mcmc)
        gelman.diag(tmp_chain_mcmc)[[1]][1,]
    })
    colnames(gelman_res_j) <- c("group","j","gelman_point","gelman_upper")
    

    data.table::setkey(n_inf_chain, "j", "group")
    n_inf_chain_summaries <- n_inf_chain[, list(
        mean = mean(as.double(total_infs)), median = median(as.double(total_infs)),
        lower_quantile = quantile(as.double(total_infs), c(0.025)),
        upper_quantile = quantile(as.double(total_infs), c(0.975)),
        effective_size = tryCatch({
            coda::effectiveSize(total_infs)
        }, error = function(e) {
            0
        })
    ),
    by = key(n_inf_chain)
    ]
    n_inf_chain_summaries <- merge(n_inf_chain_summaries, gelman_res_j, by=c("j","group"))
 
    n_inf_chain_summaries_cumu <- n_inf_chain[, list(
        mean = mean(as.double(cumu_infs)), median = median(as.double(cumu_infs)),
        lower_quantile = quantile(as.double(cumu_infs), c(0.025)),
        upper_quantile = quantile(as.double(cumu_infs), c(0.975)),
        effective_size = tryCatch({
            coda::effectiveSize(cumu_infs)
        }, error = function(e) {
            0
        })
    ),
    by = key(n_inf_chain)
    ]
    message("Done\n")
    if (!is.null(known_ar)) {
        n_inf_chain_summaries <- merge(n_inf_chain_summaries, known_ar, by = c("j","group"))
        n_inf_chain_summaries$correct <- (n_inf_chain_summaries$AR >=
                                          n_inf_chain_summaries$lower_quantile) & (n_inf_chain_summaries$AR <=
                                                                                   n_inf_chain_summaries$upper_quantile)
    }
    message("Calculating by individual summaries...\n")
    data.table::setkey(inf_chain, "i", "sampno", "chain_no")
    n_inf_chain_i <- inf_chain[, list(total_infs = sum(x)), by = key(inf_chain)]

    data.table::setkey(n_inf_chain_i, "i")
    n_inf_chain_i_summaries <- n_inf_chain_i[, list(
        mean = mean(total_infs),
        median = as.integer(median(total_infs)),
        lower_quantile = quantile(total_infs, 0.025),
        upper_quantile = quantile(total_infs, 0.975),
        effective_size = tryCatch({
            coda::effectiveSize(total_infs)
        }, error = function(e) {
            0
        })
    ),
    by = key(n_inf_chain_i)
    ]

    if(solve_cumulative){
        data.table::setkey(inf_chain,"i", "sampno", "chain_no")
        n_inf_chain_i_cumu <- inf_chain[, cumu_infs := cumsum(x), by = key(inf_chain)]

        data.table::setkey(n_inf_chain_i_cumu, "i","j")
        n_inf_chain_i_summaries_cumu <- n_inf_chain_i_cumu[, list(
            mean = mean(cumu_infs),
            median = as.integer(median(cumu_infs)),
            lower_quantile = quantile(cumu_infs, 0.025),
            upper_quantile = quantile(cumu_infs, 0.975),
            effective_size = tryCatch({
                coda::effectiveSize(cumu_infs)
            }, error = function(e) {
                0
            })
        ),
        by = key(n_inf_chain_i_cumu)
        ]
    } else {
        n_inf_chain_i_summaries_cumu <- NULL
    }
    message("Done\n")
    if (!is.null(known_infection_history)) {
        true_n_infs <- rowSums(known_infection_history)
        true_n_infs <- data.frame(i = 1:length(true_n_infs), true_infs = true_n_infs)
        n_inf_chain_i_summaries <- merge(n_inf_chain_i_summaries, true_n_infs, by = "i")
        n_inf_chain_i_summaries$correct <- (n_inf_chain_i_summaries$true_inf >=
                                            n_inf_chain_i_summaries$lower_quantile) & (n_inf_chain_i_summaries$true_inf <=
                                                                                       n_inf_chain_i_summaries$upper_quantile)
    }

    return(list(
        "by_year" = n_inf_chain_summaries, "by_indiv" = n_inf_chain_i_summaries,
        "by_year_cumu" = n_inf_chain_summaries_cumu, "by_indiv_cumu" = n_inf_chain_i_summaries_cumu
    ))
}

#' Plot historical attack rates monthly
#'
#' Plots inferred historical attack rates from the MCMC output on infection histories for monthly. The main difference compared to the normal attack rate plot is that pointrange plots don't make as much sense at a very fine time resolution.
#' @inheritParams plot_attack_rates
#' @param ymax Numeric. the maximum y value to put on the axis. Default = 1.
#' @param buckets Integer. How many buckets of time is each year split into? ie. 12 for monthly data, 4 for quarterly etc. Default = 1.
#' @param cumulative if TRUE, plots the cumulative attack rate
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates_monthly <- function(infection_histories, titre_dat, strain_isolation_times,
                                      n_alive = NULL, ymax = 1, buckets = 1,
                                      pad_chain = TRUE, true_ar = NULL, by_group = FALSE, group_subset = NULL,
                                      cumulative = FALSE,add_box=FALSE) {
    if (is.null(infection_histories$chain_no)) {
        infection_histories$chain_no <- 1
    }
    if (is.null(infection_histories$group)) {
        infection_histories$group <- 1
    }
    ## Some year/sample combinations might have no infections there.
    ## Need to make sure that these get considered
    if (pad_chain) infection_histories <- pad_inf_chain(infection_histories)
    ## Subset of groups to plot
    if (is.null(group_subset)) {
        group_subset <- unique(titre_dat$group)
    }
    if (!by_group) titre_dat$group <- 1
    ## Find inferred total number of infections from the MCMC output
    ## Scale by number of individuals that were alive in each epoch
    ## and generate quantiles
    months <- 1:buckets
    years <- strain_isolation_times
    years <- range(floor(years / buckets))
    years <- years[1]:years[2]
    labels <- c(sapply(years, function(x) paste0(months, "/", x)))
    labels1 <- labels[1:length(strain_isolation_times)]
    labels1 <- labels1[seq(1, length(labels1), by = 1)]
    year_break <- strain_isolation_times[seq(1, length(strain_isolation_times), by = 1)]

    if (is.null(n_alive)) {
        n_alive <- as.data.frame(get_n_alive_group(titre_dat, strain_isolation_times))
        n_alive$group <- 1:nrow(n_alive)
    }
    n_alive_tmp <- reshape2::melt(n_alive, id.vars = "group")
    n_alive_tmp$variable <- as.numeric(n_alive_tmp$variable)
    colnames(n_alive_tmp) <- c("group", "j", "n_alive")

    colnames(infection_histories)[1] <- "individual"
    infection_histories <- merge(infection_histories, data.table(unique(titre_dat[, c("individual", "group")])), by = c("individual","group"))
    ## Sum infections per year for each MCMC sample
    data.table::setkey(infection_histories, "sampno", "j", "group", "chain_no")
    tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]

    if (cumulative & !pad_chain) message("Error - cannot calculate cumulative incidence without pad_chain = TRUE\n")

    tmp <- merge(tmp, data.table(n_alive_tmp), by = c("group", "j"))
    tmp$time <- strain_isolation_times[tmp$j]
    tmp$V1 <- tmp$V1 / tmp$n_alive
    tmp[is.nan(tmp$V1), "V1"] <- 0
    if (cumulative && pad_chain) {
        data.table::setkey(tmp, "sampno", "group", "chain_no")
        tmp[, V1 := cumsum(V1), by = key(tmp)]
    }
    quantiles <- ddply(tmp, .(time, group), function(x) quantile(x$V1, c(0.025,0.25, 0.5,0.75, 0.975)))
    colnames(quantiles) <- c("time", "group", "lower","lower2", "median", "upper2","upper")
    
    

    p <- ggplot(quantiles[quantiles$group %in% group_subset, ])
    if(add_box){
        x_box_min <- min(titre_dat$samples)
        p <- p +geom_rect(xmin=x_box_min, xmax=max(quantiles$time),ymin=-1,ymax=2,fill="gray90",alpha=0.1)
    }
    p <- p +
        geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +
        geom_ribbon(aes(x = time, ymin = lower2, ymax = upper2), fill = "red", alpha = 0.4) +
        geom_line(aes(x = time, y = median), col = "red")
    if (!is.null(true_ar)) {
        p <- p +
            geom_line(
                data = true_ar[true_ar$group %in% group_subset, ], aes(x = j, y = AR),
                col = "purple", size = 0.5
            )
    }
    p <- p +
        ## geom_point(aes(x = year, y = median), col = "purple", size = 0.5) +
        facet_wrap(~group, ncol = 2) +
        scale_y_continuous(expand = c(0, 0)) +
                                        # #scale_x_continuous(expand = c(0, 0), breaks = year_break, labels = labels1) +
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
#' @param strain_isolation_times vector of the epochs of potential circulation
#' @param n_alive vector with the number of people alive in each year of circulation. Can be left as NULL, and ages will be used to infer this
#' @param resolution divides strain isolation times by this number for x axis labels
#' @param pointsize Numeric - how big should each point be?
#' @param fatten Numeric - fatten parameter for ggplot pointrange
#' @param pad_chain if TRUE, fills the infection history data table with entries for non-infection events (ie. 0s). Can be switched to FALSE for speed to get a rough idea of what the attack rates look like.
#' @param prior_pars if not NULL, a list of parameters for the attack rate prior, giving the assumed prior_version along with alpha and beta
#' @param plot_den if TRUE, produces a violin plot of attack rates rather than pointrange
#' @param true_ar data frame of true attack rates, with first column `year` equal to `strain_isolation_times`, and second column `AR` giving the attack rate. Column names: group, j, AR
#' @param by_group if TRUE, facets the plot by group ID
#' @param group_subset if not NULL, plots only this subset of groups eg. 1:5
#' @param plot_residuals if TRUE, plots the residuals between inferred and true attack rate
#' @param colour_by_taken if TRUE, then colours the attack rates by whether or not titres against the circulating virus at that time were measured
#' @param by_val frequency of x-axis labels
#' @return a ggplot2 object with the inferred attack rates for each potential epoch of circulation
#' @export
plot_attack_rates <- function(infection_histories, titre_dat, strain_isolation_times, 
                              n_alive = NULL,
                              pointsize = 1, fatten = 1,
                              pad_chain = TRUE, prior_pars = NULL,
                              plot_den = FALSE,
                              true_ar = NULL, by_group = FALSE,
                              group_subset = NULL, plot_residuals = FALSE,
                              colour_by_taken = TRUE, by_val = 5) {
    ## Some year/sample combinations might have no infections there.
    ## Need to make sure that these get considered
    if (is.null(infection_histories$chain_no)) {
        infection_histories$chain_no <- 1
    }

    if (pad_chain) infection_histories <- pad_inf_chain(infection_histories)

    ## Subset of groups to plot
    if (is.null(group_subset)) {
        group_subset <- unique(titre_dat$group)
    }
    if (!by_group) {
        titre_dat$group <- 1
        infection_histories$group <- 1
    }

    ## Find inferred total number of infections from the MCMC output
    ## Scale by number of individuals that were alive in each epoch
    ## and generate quantiles
    if (is.null(n_alive)) {
        n_alive <- get_n_alive_group(titre_dat, strain_isolation_times)
    }
    n_alive <- as.data.frame(n_alive)
    n_alive$group <- 1:nrow(n_alive)
    
    n_groups <- length(unique(titre_dat$group))
    n_alive_tot <- get_n_alive(titre_dat, strain_isolation_times)
    colnames(infection_histories)[1] <- "individual"
    infection_histories <- merge(infection_histories, data.table(unique(titre_dat[, c("individual", "group")])), by = c("individual","group"))
    years <- c(strain_isolation_times, max(strain_isolation_times) + 3)
    data.table::setkey(infection_histories, "sampno", "j", "chain_no", "group")
    tmp <- infection_histories[, list(V1 = sum(x)), by = key(infection_histories)]
    tmp$taken <- years[tmp$j] %in% unique(titre_dat$samples)
    tmp$taken <- ifelse(tmp$taken, "Yes", "No")
    prior_dens <- NULL
    n_alive1 <- n_alive
    if (!is.null(prior_pars)) {
        n_alive$Prior <- 1
        prior_ver <- prior_pars[["prior_version"]]
        alpha1 <- prior_pars[["alpha"]]
        beta1 <- prior_pars[["beta"]]
        if (prior_ver == 3) {
            prior_dens <- rbinom(10000, size = max(n_alive_tot), p = alpha1 / (alpha1 + beta1)) / max(n_alive_tot)
        } else {
            prior_dens <- rbeta(10000, alpha1, beta1)
        }
        prior_dens <- data.frame(
            sampno = 1:length(prior_dens), j = max(tmp$j) + 1,
            chain_no = 1, V1 = prior_dens, taken = "Prior", group = 1
        )
        prior_dens_all <- NULL
        for (i in 1:nrow(n_alive)) {
            prior_dens$group <- i
            prior_dens_all <- rbind(prior_dens_all, prior_dens)
        }
        tmp <- rbind(tmp, prior_dens_all)
    }
    n_alive_tmp <- reshape2::melt(n_alive, id.vars = "group")
    n_alive_tmp$variable <- as.numeric(n_alive_tmp$variable)
    colnames(n_alive_tmp) <- c("group", "j", "n_alive")
    tmp <- merge(tmp, data.table(n_alive_tmp), by = c("group", "j"))
    tmp$V1 <- tmp$V1 / tmp$n_alive

    min_year <- min(strain_isolation_times)
    max_year <- max(strain_isolation_times)
    year_breaks <- c(min_year, seq(5 * round(min_year / 5), max_year, by = 5))
    year_labels <- c(min_year, seq(5 * round(min_year / 5), max_year, by = 5))

    if (!plot_den) {
        quantiles <- ddply(tmp, .(j, group), function(x) quantile(x$V1, c(0.025, 0.5, 0.975)))
        colnames(quantiles) <- c("j", "group", "lower", "median", "upper")
                                        # quantiles[c("lower", "median", "upper")] <- quantiles[c("lower", "median", "upper")]# / n_alive1[quantiles$j]
        quantiles$j <- years[quantiles$j]
        quantiles$taken <- quantiles$j %in% unique(titre_dat$samples)
        quantiles$taken <- ifelse(quantiles$taken, "Yes", "No")

        quantiles$tested <- quantiles$j %in% unique(titre_dat$virus)
        quantiles$tested <- ifelse(quantiles$tested, "Yes", "No")

        min_year <- min(strain_isolation_times)
        max_year <- max(strain_isolation_times)
        year_breaks <- c(min_year, seq(by_val * round(min_year / by_val), max_year, by = by_val))
        year_labels <- c(min_year, seq(by_val * round(min_year / by_val), max_year, by = by_val))

        if (!is.null(prior_dens)) {
            quantiles[quantiles$j == max(years), "taken"] <- "Prior"
            year_breaks <- c(year_breaks, max_year + 3)
            year_labels <- c(year_labels, "Prior")
        }
        colnames(quantiles)[which(colnames(quantiles) == "taken")] <- "Sample taken"
        colnames(quantiles)[which(colnames(quantiles) == "tested")]  <- "Virus tested"

        p <- ggplot(quantiles[quantiles$group %in% group_subset, ]) +
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
                         col = `Virus tested`, shape = `Virus tested`
                     ),
                     size = pointsize,
                     fatten = fatten
                     )
        }
    } else {
        tmp$j <- years[tmp$j]

        if (plot_residuals) {
            true_ar <- true_ar[, c("group", "j", "AR")]
            if (!is.null(prior_pars)) {
                true_ar <- rbind(true_ar, data.frame(
                                              group = 1:n_groups,
                                              j = max(strain_isolation_times) + 3,
                                              AR = median(prior_dens$V1)
                                          ))
            }
            tmp <- merge(tmp, true_ar, by = c("group", "j"))
            tmp$V1 <- tmp$V1 - tmp$AR
        }

        p <- ggplot(tmp[tmp$group %in% group_subset, ]) +
            geom_violin(aes(x = j, y = V1, fill = taken, group = j),
                        draw_quantiles = c(0.5), scale = "width",
                        adjust=2
                        )
    }
    if (!is.null(true_ar) & !plot_residuals) {
        p <- p +
            geom_point(
                data = true_ar[true_ar$group %in% group_subset, ], aes(x = j, y = AR),
                col = "purple", size = 0.5
            )
    }
    if (!plot_residuals) {
        p <- p +
            scale_y_continuous(expand = c(0, 0)) +
            coord_cartesian(ylim=c(0,1))
    }

    if (by_group) {
        p <- p + facet_wrap(~group, ncol = 2)
    }
    p <- p +
        scale_x_continuous(breaks = year_breaks, labels = year_labels) +
        theme_classic() +
        theme(legend.position = "none") +
        ylab("Estimated attack rate") +
        xlab("Year")
    if (plot_residuals) {
        p <- p +
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_y_continuous(limits = c(-1, 1), expand = c(0, 0))
        p_res <- ggplot(tmp) + geom_density(aes(x = V1), fill = "grey40") +
            facet_wrap(~group) +
            geom_vline(xintercept = 0, linetype = "dashed", colour = "blue") +
            scale_x_continuous(limits = c(-1, 1)) +
            theme_classic() + xlab("Estimated AR - true AR") + ylab("Density")
        return(list(p, p_res))
    }

    return(p)
}

#' Pad infection history chain
#'
#' Given that the MCMC sampler only stores present infections (ie. there are no entries for 0s from the infection history matrix), for some summaries we need to add these 0s back in to avoid bias.
#' @param inf_chain the data table with infection history samples from \code{\link{run_MCMC}}
#' @return the same inf_chain that was passed in, but with 0s for missing i/j/sampno combinations
#' @export
pad_inf_chain <- function(inf_chain) {
    if (is.null(inf_chain$chain_no)) {
        inf_chain$chain_no <- 1
    }
    if (is.null(inf_chain$group)) {
        inf_chain$group <- 1
    }
    
    is <- unique(inf_chain$i)
    js <- unique(inf_chain$j)

    sampnos <- unique(inf_chain$sampno)
    chain_nos <- unique(inf_chain$chain_no)
    groups <- unique(inf_chain$group)
    expanded_values <- data.table::CJ(
                                       i = is,
                                       j = js,
                                       sampno = sampnos,
                                       chain_no = chain_nos,
                                       group = groups
                                   )
    diff_infs <- fsetdiff(expanded_values, inf_chain[, c("i", "j", "sampno", "chain_no", "group")])
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
#' @param pad_chain if TRUE, pads the infection history MCMC chain to have entries for non-infection events
#' @return a list of two ggplot objects - the MCMC trace and MCMC densities
#' @seealso \code{\link{plot_infection_history_chains_indiv}}
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_inf_chain)
#' data(example_titre_dat)
#' data(example_antigenic_map)
#' times <- example_antigenic_map$inf_times
#' n_alive_group <- get_n_alive_group(example_titre_dat, strain_isolation_times,melt_dat = TRUE)
#' n_alive_group$j <- strain_isolation_times[n_alive_group$j]
#' plot_infection_history_chains_time(example_inf_chain, 0, sample(1:length(times),10),n_alive,FALSE)
#' }
#' @export
plot_infection_history_chains_time <- function(inf_chain, burnin = 0, years = NULL,
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
#' plot_number_infections(example_inf_chain)
#' }
#' @export
plot_number_infections <- function(inf_chain, pad_chain = TRUE) {
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
                        size = 0.1, shape = 21, fatten = 0.1
                        ) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab("Individual (ordered)") +
        ylab("Estimated number of infections") +
        theme_bw()
    return(p)
}


#' Useful plot for looking at simulated data
#'
#' Plots measured titres and known infection histories for all individuals, facetted by sample time (multi-strain panel) or virus variable (longitidunal single strain)
#' @param titre_dat the data frame of titre data
#' @param infection_histories the infection history matrix
#' @param strain_isolation_times the vector of times at which individuals could be infected
#' @param n_indivs how many individuals to plot
#' @param start_inf if not NULL, plots the infection history matrix used as the starting point in the MCMC chain
#' @param study_design default "multi-strain" facets by sample time. "single-strain" gives sample time on the x-axis and colours by virus
#' @return a ggplot object
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_titre_dat)
#' data(example_inf_hist)
#' data(example_antigenic_map)
#' strain_isolation_times <- example_antigenic_map$inf_times
#' plot_data(example_titre_dat, example_inf_hist, strain_isolation_times, 5)
#' }
#' @export
plot_data <- function(titre_dat, infection_histories, 
                      strain_isolation_times, 
                      n_indivs, start_inf = NULL,
                      study_design="multi-strain"){
    indivs <- unique(titre_dat$individual)
    infection_history <- as.data.frame(cbind(indivs, infection_histories))
    colnames(infection_history) <- c("individual", strain_isolation_times)
    melted_inf_hist <- reshape2::melt(infection_history, id.vars = "individual")
    melted_inf_hist$variable <- as.numeric(as.character(melted_inf_hist$variable))
    melted_inf_hist <- melted_inf_hist[melted_inf_hist$value > 0, ]
    tmp <- unique(titre_dat[, c("individual", "samples")])
    melted_inf_hist <- merge(melted_inf_hist, tmp)
    melted_inf_hist <- melted_inf_hist[melted_inf_hist$variable <= melted_inf_hist$samples, ]
    samps <- sample(unique(titre_dat$individual), n_indivs)

    if (study_design == "multi-strain") {
        p1 <- ggplot(titre_dat[titre_dat$individual %in% samps, ]) +
            geom_point(aes(x = as.integer(virus), y = titre)) +
            geom_vline(data = melted_inf_hist[melted_inf_hist$individual %in% samps, ], aes(xintercept = variable), col = "red", linetype = "dashed") +
            theme_bw() +
            xlab("Strain") +
            facet_grid(individual ~ samples)
    } else {
        p1 <- ggplot(titre_dat[titre_dat$individual %in% samps, ]) +
            geom_point(aes(x = samples, y = titre, col=virus)) +
            geom_vline(data = melted_inf_hist[melted_inf_hist$individual %in% samps, ], 
                       aes(xintercept = variable), col = "red", linetype = "dashed") +
            theme_bw() +
            xlab("Strain circulation time") +
            facet_wrap(~individual)
    }
    
    if (!is.null(start_inf)) {
        start_inf_hist <- as.data.frame(cbind(indivs, start_inf))
        colnames(start_inf_hist) <- c("individual", strain_isolation_times)
        melted_start_hist <- reshape2::melt(start_inf_hist, id.vars = "individual")
        melted_start_hist$variable <- as.numeric(as.character(melted_start_hist$variable))
        melted_start_hist <- melted_start_hist[melted_start_hist$value > 0, ]
        p1 <- p1 + geom_vline(data = melted_start_hist[melted_start_hist$individual %in% samps, ], aes(xintercept = variable), col = "blue", linetype = "dashed")
    }
    p1 <- p1 + ylab("log titre")
    return(p1)
}

#' Plot cumulative and per time posterior infection probability densities
#'
#' For each individual requested, plots the median and 95% quantiles on a) the cumulative number of infections over a lifetime and b) the posterior probability that an infection occured in a given time point
#' @param inf_chain the infection history chain
#' @param burnin only plot samples where sampno > burnin
#' @param indivs vector of individual ids to plot
#' @param real_inf_hist if not NULL, adds lines to the plots showing the known true infection times
#' @param start_inf if not NULL, adds lines to show where the MCMC chain started
#' @param strain_isolation_times vector of times at which individuals could have been infected
#' @param nsamp how many samples from the MCMC chain to take?
#' @param ages if not NULL, adds lines to show when an individual was born
#' @param number_col how many columns to use for the cumulative infection history plot
#' @param pad_chain if TRUE, pads the infection history MCMC chain to have entries for non-infection events
#' @param subset_years if not NULL, pass a vector of indices to only take a subset of indices from strain_isolation_times
#' @param return_data if TRUE, returns the infection history posterior densities used to generate the plots
#' @return two ggplot objects
#' @family infection_history_plots
#' @examples
#' \dontrun{
#' data(example_inf_chain)
#' data(example_antigenic_map)
#' data(example_inf_hist)
#' data(example_titre_dat)
#' 
#' ages <- unique(example_titre_dat[,c("individual","DOB")])
#' times <- example_antigenic_map$inf_times
#' indivs <- 1:10
#' generate_cumulative_inf_plots(example_inf_chain, 0, indivs, example_inf_hist, NULL, times,
#'                               ages=ages, number_col=2,pad_chain=FALSE, return_data=TRUE)
#' }
#' @export
generate_cumulative_inf_plots <- function(inf_chain, burnin = 0, indivs, real_inf_hist = NULL, start_inf = NULL,
                                          strain_isolation_times, nsamp = 100, ages = NULL, number_col = 1,
                                          pad_chain = TRUE, subset_years = NULL, return_data = FALSE) {
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
    if (!is.null(subset_years)) inf_chain1 <- inf_chain1[inf_chain1$j %in% subset_years, ]
    data.table::setkey(inf_chain1, "i", "j", "chain_no")
                                        # max_sampno <- length(unique(inf_chain1$sampno))
    max_sampno <- nsamp

    ## Number of samples with a 1 divided by total samples
    densities <- inf_chain1[, list(V1 = sum(x) / max_sampno), by = key(inf_chain1)]

    ## Convert to real time points
    densities$j <- as.numeric(strain_isolation_times[densities$j])
    densities$i <- as.numeric(densities$i)

    ## If someone wasn't infected in a given year at all, then need a 0
    strain_isolation_times1 <- strain_isolation_times
    if (!is.null(subset_years)) strain_isolation_times1 <- strain_isolation_times[subset_years]
    all_combos <- data.table(expand.grid(i = indivs, j = strain_isolation_times1, chain_no = unique(inf_chain$chain_no)))
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

        colnames(infection_history1) <- c(strain_isolation_times, "i")
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
        xlab("Year") +
        ylab("Density") +
        theme(legend.position = "bottom") +
        theme_bw()

    ## Generate lower, upper and median cumulative infection histories from the
    ## MCMC chain
    tmp_inf_chain <- inf_chain[inf_chain$i %in% indivs, ]
    hist_profiles <- ddply(tmp_inf_chain, .(i, sampno, chain_no), function(x) {
        empty <- numeric(length(strain_isolation_times))
        empty[x[x$x == 1, "j"]] <- 1
        cumsum(empty)
    })

    hist_profiles <- hist_profiles[, colnames(hist_profiles) != "sampno"]

    colnames(hist_profiles) <- c("i", "chain_no", strain_isolation_times)
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
        age_mask <- create_age_mask(tmp_age[, 2], strain_isolation_times)
        age_dat <- data.frame(j = age_mask, individual = indivs[order(indivs)])
        p1 <- p1 + geom_vline(data = age_dat, aes(xintercept = strain_isolation_times[j]), col = "purple", linetype = "dashed")
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
#' @param titre_dat the data frame of titre data, including labels for individuals and time sample was taken
#' @return a ggplot2 object
#' @family theta_plots
#' @examples
#' \dontrun{
#' data(example_titre_dat)
#' plot_samples_distances(example_titre_dat)
#' }
#' @export
plot_samples_distances <- function(titre_dat) {
    samples <- unique(titre_dat[, c("individual", "samples")])
    distances <- ddply(samples, ~individual, function(x) {
        if (nrow(x) < 2) {
            y <- 0
        } else {
            y <- diff(x$samples)
        }
        y
    })
    ggplot(distances) + geom_histogram(aes(x = V1), binwidth = 1) + theme_bw() + xlab("Time points between samples")
}
