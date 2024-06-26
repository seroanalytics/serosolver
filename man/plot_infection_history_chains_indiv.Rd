% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_infection_histories.R
\name{plot_infection_history_chains_indiv}
\alias{plot_infection_history_chains_indiv}
\title{Plot MCMC trace for infections per individual}
\usage{
plot_infection_history_chains_indiv(
  inf_chain,
  burnin = 0,
  indivs = NULL,
  pad_chain = TRUE
)
}
\arguments{
\item{inf_chain}{the data table with infection history samples from \code{\link{serosolver}}}

\item{burnin}{optionally remove all samp_no < burnin from the chain}

\item{indivs}{vector of integers, if not NULL, only plots a subset of individuals (where 1 is the first individual)}

\item{pad_chain}{if TRUE, pads the infection history MCMC chain to have entries for non-infection events}
}
\value{
a list of two ggplot objects - the MCMC trace and MCMC densities
}
\description{
Plot MCMC trace for infections per individual
}
\examples{
\dontrun{
data(example_inf_chain)
plot_infection_history_chains_indiv(example_inf_chain, 0, 1:10, FALSE)
}
}
\seealso{
\code{\link{plot_infection_history_chains_indiv}}

Other infection_history_plots: 
\code{\link{calculate_infection_history_statistics}()},
\code{\link{plot_antibody_data}()},
\code{\link{plot_antibody_predictions}()},
\code{\link{plot_cumulative_infection_histories}()},
\code{\link{plot_estimated_antibody_model}()},
\code{\link{plot_individual_number_infections}()},
\code{\link{plot_infection_history_chains_time}()},
\code{\link{plot_model_fits}()},
\code{\link{plot_posteriors_infhist}()},
\code{\link{plot_total_number_infections}()}
}
\concept{infection_history_plots}
