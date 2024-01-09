
#' Calculate FOI log probability
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @param phis a vector of FOIs
#' @param infection_history the matrix of infection histories
#' @param age_mask the age mask vector as returned by \code{\link{create_age_mask}}
#' @return a single log probability
#' @family priors
#' @export
calc_phi_probs <- function(phis, infection_history, age_mask, sample_mask) {
  lik <- 0
  for (i in 1:ncol(infection_history)) {
    use_indivs <- intersect(which(age_mask <= i), which(sample_mask >= i))
    lik <- lik + sum(log((phis[i]^infection_history[use_indivs, i] *
                            (1 - phis[i])^(1 - infection_history[use_indivs, i]))))
  }
  lik
}

#' Calculate FOI log probability vector
#'
#' Given a vector of FOIs for all circulating years, a matrix of infection histories and the vector specifying if individuals were alive or not, returns the log probability of the FOIs given the infection histories.
#' @inheritParams calc_phi_probs
#' @return a vector of log probabilities for each individual
#' @family priors
#' @export
calc_phi_probs_indiv <- function(phis, infection_history, age_mask, sample_mask) {
  lik <- numeric(nrow(infection_history))
  for (i in 1:ncol(infection_history)) {
    to_add <- log(((phis[i]^infection_history[, i]) *
                     (1 - phis[i])^(1 - infection_history[, i]))) *
      as.numeric(age_mask <= i) *
      as.numeric(sample_mask >= i)
    lik <- lik + to_add
  }
  lik
}

#' Calculate FOI from spline - INACTIVE
#'
#' Version of FOI prior that calculates a seasonal spline such that the FOI in a given time point (month) comes from this spline term
#' @inheritParams calc_phi_probs
#' @param foi vector of force of infections per year
#' @param knots vector of knots for the spline
#' @param theta vector of theta parameters for spline
#' @return a vector of log probabilities for each individual
#' @family priors
calc_phi_probs_spline <- function(foi, knots, theta, infection_history, age_mask) {
  phis <- generate_phis(foi, knots, theta, length(foi), 12)
  lik <- numeric(nrow(infection_history))
  for (i in 1:ncol(infection_history)) {
    lik <- lik + infection_history[, i] * log(phis[i]) + (1 - infection_history[, i]) * log(phis[i]) + log(as.numeric(age_mask <= i)) + log(as.numeric(sample_mask >= i))
  }
  lik
}

#' Generate FOI phis - INACTIVE
#'
#' Calculates the seasonal spline for \code{\link{calc_phi_probs_spline}}
#' @inheritParams calc_phi_probs_spline
#' @param n_years number of years to calculate
#' @param buckets number of buckets per year (12 for monthly, 1 for annual)
#' @param degree degree of the spline
#' @return a vector of FOIs for each time point
generate_phis <- function(foi, knots, theta, n_years, buckets, degree = 2) {
  x <- seq(0, buckets - 1, by = 1) / buckets
  n_knots <- length(knots) + degree + 1
  all_dat <- NULL
  index <- 1
  tmp <- gen_spline_y(x, knots, degree, theta)
  tmp <- tmp / sum(tmp)
  all_dat <- numeric(length(tmp) * length(foi))
  for (i in 1:n_years) {
    all_dat[index:(index + length(tmp) - 1)] <- tmp * foi[i]
    index <- index + length(tmp)
  }
  all_dat
}

#' Generates a spline for \code{\link{generate_phis}} - INACTIVE
gen_spline_y <- function(x, knots, degree, theta, intercept = TRUE) {
  basis <- bs(
    x = x, knots = knots, degree = degree,
    Boundary.knots = c(0, 1), intercept = intercept
  )
  
  y.spline <- basis %*% theta
  return(as.vector(y.spline))
}

calc_phi_loc_probs_indiv <- function(phis, group_probs, infection_history, age_mask, sample_mask, group_indices) {
  lik <- numeric(nrow(infection_history))
  max_group_p <- max(group_probs)
  group_probs <- group_probs / max_group_p
  for (i in 1:ncol(infection_history)) {
    tmp <- (infection_history[, i] * log(phis[i] * group_probs[group_indices]) +
              (1 - infection_history[, i]) * log(1 - phis[i] * group_probs[group_indices])) *
      as.numeric(age_mask <= i) * as.numeric(sample_mask >= i)
    lik <- lik + tmp
  }
  lik
}