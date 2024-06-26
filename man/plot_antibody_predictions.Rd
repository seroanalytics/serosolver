% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_antibody_model.R
\name{plot_antibody_predictions}
\alias{plot_antibody_predictions}
\title{Plots model predicted titers against observations}
\usage{
plot_antibody_predictions(
  chain,
  infection_histories,
  antibody_data = NULL,
  demographics = NULL,
  par_tab = NULL,
  antigenic_map = NULL,
  possible_exposure_times = NULL,
  nsamp = 1000,
  measurement_bias = NULL,
  data_type = 1,
  start_level = "none",
  settings = NULL
)
}
\arguments{
\item{chain}{the full MCMC chain to generate antibody level trajectories from}

\item{infection_histories}{the MCMC chain for infection histories}

\item{antibody_data}{the data frame of antibody level data}

\item{par_tab}{the table controlling the parameters in the MCMC chain}

\item{antigenic_map}{(optional) a data frame of antigenic x and y coordinates. Must have column names: x_coord; y_coord; inf_times. See \code{\link{example_antigenic_map}}}

\item{possible_exposure_times}{(optional) if no antigenic map is specified, this argument gives the vector of times at which individuals can be infected}

\item{nsamp}{number of samples to take from posterior}

\item{measurement_bias}{default NULL, optional data frame giving the index of `rho` that each biomarker_id and biomarker_group which uses the measurement shift from from. eg. if there's 6 circulation years and 3 strain clusters}

\item{data_type}{integer, currently accepting 1 or 2. Set to 1 for discretized, bounded data, or 2 for continuous, bounded data.}

\item{settings}{if not NULL, list of serosolver settings as returned from the main serosolver function}
}
\value{
a list with: 
\itemize{
\item a data frame with all posterior estimates for each observation; 
\item the proportion of observations captured by the 95% prediction intervals; 
\item a histogram comparing posterior median estimates to the observed data (note, this can be misleading for continuous data due to the zero-inflated observation model); 
\item a histogram comparing random posterior draws to the observed data (can be more reliable than posterior medians);
\item comparison of observations and all posterior medians and 95% prediction intervals
}
}
\description{
Plots model predicted titers against observations
}
\seealso{
Other infection_history_plots: 
\code{\link{calculate_infection_history_statistics}()},
\code{\link{plot_antibody_data}()},
\code{\link{plot_cumulative_infection_histories}()},
\code{\link{plot_estimated_antibody_model}()},
\code{\link{plot_individual_number_infections}()},
\code{\link{plot_infection_history_chains_indiv}()},
\code{\link{plot_infection_history_chains_time}()},
\code{\link{plot_model_fits}()},
\code{\link{plot_posteriors_infhist}()},
\code{\link{plot_total_number_infections}()}
}
\concept{infection_history_plots}
