% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_infection_histories.R
\name{plot_individual_number_infections}
\alias{plot_individual_number_infections}
\title{Plot point range number infections per individual}
\usage{
plot_individual_number_infections(inf_chain, pad_chain = TRUE)
}
\arguments{
\item{inf_chain}{the data table with infection history samples from \code{\link{serosolver}}}

\item{pad_chain}{if TRUE, pads the infection history MCMC chain to have entries for non-infection events}
}
\value{
a ggplot object
}
\description{
Plot point range number infections per individual
}
\examples{
\dontrun{
data(example_inf_chain)
plot_individual_number_infections(example_inf_chain)
}
}
\seealso{
Other infection_history_plots: 
\code{\link{calculate_infection_history_statistics}()},
\code{\link{plot_antibody_data}()},
\code{\link{plot_antibody_predictions}()},
\code{\link{plot_cumulative_infection_histories}()},
\code{\link{plot_estimated_antibody_model}()},
\code{\link{plot_infection_history_chains_indiv}()},
\code{\link{plot_infection_history_chains_time}()},
\code{\link{plot_model_fits}()},
\code{\link{plot_posteriors_infhist}()},
\code{\link{plot_total_number_infections}()}
}
\concept{infection_history_plots}
