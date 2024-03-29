% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{inf_mat_prior_cpp}
\alias{inf_mat_prior_cpp}
\title{Marginal prior probability (p(Z)) of a particular infection history matrix single prior
 Prior is independent contribution from each year}
\usage{
inf_mat_prior_cpp(infection_history, n_alive, shape1, shape2)
}
\arguments{
\item{infection_history}{IntegerMatrix, the infection history matrix}

\item{n_alive}{IntegerVector, vector giving the number of individuals alive in each time unit}

\item{shape1}{double, shape1 (alpha) parameter for beta distribution prior}

\item{shape2}{double, shape2 (beta) parameter for beta distribution prior}
}
\value{
a single prior probability
}
\description{
Marginal prior probability (p(Z)) of a particular infection history matrix single prior
 Prior is independent contribution from each year
}
\seealso{
Other inf_mat_prior: 
\code{\link{inf_mat_prior_cpp_vector}()},
\code{\link{inf_mat_prior_group_cpp_vector}()},
\code{\link{inf_mat_prior_group_cpp}()},
\code{\link{inf_mat_prior_total_group_cpp}()}
}
\concept{inf_mat_prior}
