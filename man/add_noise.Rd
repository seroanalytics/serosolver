% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{add_noise}
\alias{add_noise}
\title{Add noise}
\usage{
add_noise(y, theta, measurement_bias = NULL, indices = NULL, data_type = 1)
}
\arguments{
\item{y}{the titre}

\item{theta}{a vector with max_measurement and error parameters}

\item{data_type}{integer, currently accepting 1 or 2. Set to 1 for discretized, bounded data, or 2 for continuous, bounded data. 3 is for continuous data assuming that true negatives follow a different distribution -- the vast majority return the min_measurement, but with a rate fp_rate, a random draw from a uniform distribution within the limits of detection is generated.}
}
\value{
a noisy titre
}
\description{
Adds truncated noise to antibody data
}
\examples{
\dontrun{
## ... example in simulate_individual
pars <- c("obs_sd"=1)
y <- runif(100)
noisy_y <- add_noise(y, pars)
}
}
