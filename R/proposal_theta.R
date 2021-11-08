#' MCMC proposal function
#'
#' Proposal function for MCMC random walk, taking random steps of a given size.
#' @param values a vector of the parameters to be explored
#' @param lower_bounds a vector of the low allowable bounds for the proposal
#' @param upper_bounds a vector of the upper allowable bounds for the proposal
#' @param steps a vector of step sizes for the proposal
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @family proposals
#' @export
#' @useDynLib serosolver
univ_proposal <- function(values, lower_bounds, upper_bounds, steps, index) {
  rtn <- values

  ## Commented out using guassian proposals, which does work but risks being
  ## inefficient due to stepping outside of bounds
  ## rtn[index] <- rnorm(1,values[index],steps[index])
  ## return(rtn)

  mn <- lower_bounds[index]
  mx <- upper_bounds[index]

  rtn <- values
  x <- toUnitScale(values[index], mn, mx)

  stp <- steps[index]

  rv <- runif(1)
  rv <- (rv - 0.5) * stp
  x <- x + rv

  ## Bouncing boundary condition
  if (x < 0) x <- -x
  if (x > 1) x <- 2 - x

  ## Cyclical boundary conditions
  # if (x < 0) x <- 1 + x
  # if (x > 1) x <- x - 1

  if (x < 0 | x > 1) {
    print("Stepped outside of unit scale. Something went wrong...")
  }

  rtn[index] <- fromUnitScale(x, mn, mx)
  rtn
}


#' Multivariate proposal function
#'
#' Given the current parameters and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution. Takes into account parameter covariance and ensures containment condition with beta, if cov_mat0 (the identity matrix) is specified.
#' @param values the vector of current parameter values
#' @param fixed set of flags corresponding to the parameter vector indicating which parameters are fixed
#' Takes into account parameter covariance and ensures containment condition with beta, if cov_mat0 (the identity matrix) is specified.
#' @param cov_mat the 2D covariance matrix for all of the parameters
#' @param cov_mat0 optional, usually the identity matrix for theta
#' @param use_log flag. If TRUE, propose on log scale
#' @param beta Beta as in Rosenthal and Roberts 2009
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @family proposals
#' @export
#' @useDynLib serosolver
mvr_proposal <- function(values, fixed, cov_mat, cov_mat0 = NULL, use_log = FALSE, beta = 0.05) {
  proposed <- values
  ## Sample either from a single covariance matrix or weighted average of the identity matrix and
  ## given cov matrix, if specified. On either a log or linear scale.
  if (is.null(cov_mat0)) {
    if (!use_log) {
      proposed[fixed] <- MASS::mvrnorm(n = 1, mu = proposed[fixed], Sigma = (5.6644 / length(fixed)) * cov_mat)
    } else {
      proposed[fixed] <- exp(MASS::mvrnorm(n = 1, mu = log(proposed[fixed]), Sigma = (5.6644 / length(fixed)) * cov_mat))
    }
  } else {
    if (!use_log) {
      proposed[fixed] <-
        (1 - beta) * MASS::mvrnorm(n = 1, mu = proposed[fixed], Sigma = (5.6644 / length(fixed)) * cov_mat) +
        beta * MASS::mvrnorm(n = 1, mu = proposed[fixed], Sigma = (0.01 / length(fixed)) * cov_mat0)
    } else {
      proposed[fixed] <-
        (1 - beta) * exp(MASS::mvrnorm(n = 1, mu = log(proposed[fixed]), Sigma = (5.6644 / length(fixed)) * cov_mat)) +
        beta * exp(MASS::mvrnorm(n = 1, mu = log(proposed[fixed]), Sigma = (0.01 / length(fixed)) * cov_mat0))
    }
  }
  return(proposed)
}
