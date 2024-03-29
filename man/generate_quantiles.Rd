% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{generate_quantiles}
\alias{generate_quantiles}
\title{Formatted quantiles}
\usage{
generate_quantiles(x, sig_f = 3, qs = c(0.025, 0.5, 0.975), as_text = TRUE)
}
\arguments{
\item{x}{the vector to summarise}

\item{sig_f}{how many significant figures to print}

\item{qs}{the vector of quantiles}

\item{as_text}{if TRUE, formats nicely as text rather than a vector of numbers}
}
\value{
the formatted quantiles
}
\description{
Given a vector of MCMC samples, generates and formats the desired quantile estimates
}
\examples{
data(example_theta_chain)
x <- example_theta_chain$boost_long
generate_quantiles(x)
}
