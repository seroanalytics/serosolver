% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/observation_model_likelihoods.R
\name{prob_shifts}
\alias{prob_shifts}
\title{Prior on measurement shifts}
\usage{
prob_shifts(rhos, pars)
}
\arguments{
\item{rhos}{vector of measurement shifts}

\item{pars}{vector of parameters, including rho_mean and rho_sd for the normal distribution}
}
\value{
a single prior probability
}
\description{
Assumes measurement shifts are drawn from a normal distribution (random effects) with given standard deviation and mean. Code is commented out to assume log normally distributied
}
\seealso{
Other priors: 
\code{\link{calc_phi_probs_indiv}()},
\code{\link{calc_phi_probs_spline}()},
\code{\link{calc_phi_probs}()},
\code{\link{create_prob_shifts}()},
\code{\link{find_beta_prior_mode}()},
\code{\link{find_beta_prior_with_mean_var}()},
\code{\link{find_beta_prior_with_mean}()},
\code{\link{fit_beta_prior}()},
\code{\link{fit_normal_prior}()},
\code{\link{inf_mat_prior}()},
\code{\link{infection_history_prior}()}
}
\concept{priors}
